import Foundation

class Equation : Codable {
    var T: Double
    var K: Int
    var M: Int
    var N: Int
    var R: Double
    
    var h_r: Double
    var h_phi: Double
    var tau: Double
    
    init(_T: Double, _K: Int, _M: Int, _N: Int, _R: Double, _h_r: Double, _h_phi: Double, _tau: Double){
        T = _T
        K = _K
        M = _M
        N = _N
        R = _R
        
        h_r = _h_r
        h_phi = _h_phi
        tau = _tau
    }
    
    private func Calculate_f(w: [Double]) -> Array<Array<[Double]>> {
        var f = Array(repeating: Array(repeating: Array(repeating: 0.0, count: M), count: K+2), count: N)
        
        for n in 0...N-1 {
            for k in 1...K {
                //for m in 0...M {
                f[n][k][Settings.ActuatorIndex] = w[n]
                //Вариант для сравнения с аналитическим значением
                //f[n][k][m] = Settings.Example(type: Settings.FunctionType.f, m: m, n: Double(n))
                //}
            }
        }
        
        return f
    }
    
    private func Calculate_h0() -> Array<[Double]> {
        var h0 = Array(repeating: Array(repeating: 0.0, count: M), count: K+2)
        
        for k in 1...K {
            for m in 0...M-1 {
                h0[k][m] = Settings.Example(type: Settings.FunctionType.h0, m: m)
            }
        }
        
        return h0
    }
    
    private func Calculate_h1() -> Array<[Double]> {
        var h1 = Array(repeating: Array(repeating: 0.0, count: M), count: K+2)
        
        for k in 1...K {
            for m in 0...M-1 {
                h1[k][m] = Settings.Example(type: Settings.FunctionType.h1, m: m)
            }
        }
        
        return h1
    }
    
    private func GetFirstTimeLayer(prev: Array<[Double]>, f: Array<[Double]>, h1: Array<[Double]>) -> Array<[Double]> {
        var next = Array(repeating: Array(repeating: 0.0, count: M), count: K+2)
        // next[k][M] остаётся равным нулю
        
        for k in 1...K {
            //u_km^1
            for m in 0...M-2 {
                next[k][m] =    prev[k][m] +
                    tau * h1[k][m] +
                    0.5 * pow(tau / h_r, 2) * Braces(prev: prev, f: f, k: k, m: m)
            }
        }
        
        return next
    }
    
    private func GetOtherTimeLayers(cur: Array<[Double]>, prev: Array<[Double]>, f: Array<[Double]>) -> Array<[Double]> {
        // next[k][M] остаётся равным нулю
        var next = Array(repeating: Array(repeating: 0.0, count: M), count: K+2)
        
        //u_km^n+1
        for k in 1...K {
            for m in 0...M-2 {
                next[k][m] =    2 * cur[k][m] -
                    prev[k][m] +
                    pow(tau / h_r, 2) * Braces(prev: prev, f: f, k: k, m: m)
            }
        }
        
        return next
    }
    
    private func Braces(prev: Array<[Double]>, f: Array<[Double]>, k: Int, m: Int) -> Double {
        
        let temp1 = (Double(m) + 1) * (prev[k][m+1] - prev[k][m])
        let temp2 = m == 0 ? 0 : Double(m) * (prev[k][m] - prev[k][m - 1])
        let const1 = (temp1 - temp2) / (Double(m) + 0.5)
        let const2 = (prev[k + 1][m] - 2 * prev[k][m] + prev[k - 1][m]) / pow((Double(m) + 0.5) * h_phi, 2)
        
        return  (const1 +
            const2) +
            f[k][m] * pow(h_r, 2)
    }
    
    public func Solve(w: [Double]) -> Array<Array<[Double]>> {
        var solution = Array(repeating: Array(repeating: Array(repeating: 0.0, count: M), count: K+2), count: N)
        
        let f = Calculate_f(w: w)
        
        for n in 0...N-1 {
            switch(n) {
            case 0:
                solution[n] = Calculate_h0()
                break;
            case 1:
                solution[n] = GetFirstTimeLayer(prev: solution[0], f: f[0], h1: Calculate_h1())
                break;
            default:
                solution[n] = GetOtherTimeLayers(cur: solution[n - 1], prev: solution[n - 2], f: f[n - 1])
                break;
            }
            
            shift(u: &solution[n])
        }
        
        
        //print("Решение: \(solution)")
        //let maxDiff = MathHelper.calculateMaxDiff(actual: solution, expected: MathHelper.calculateExactFunction())
        //print("exact: \(MathHelper.calculateExactFunction()[N])")
        //print("diff: \(maxDiff)")
        //let maxDiff = MathHelper.calculateMaxDiff1(actual: solution, expected: MathHelper.calculateExactFunction())
        //print("exact: \(MathHelper.calculateExactFunction()[N])")
        //print("diff: \(maxDiff)")
        return solution
    }
    
    private func shift(u: inout Array<[Double]>) {
        for m in 0...M-1 {
            u[0][m] = u[2][m]
            u[K + 1][m] = u[K - 1][m]
        }
    }
    
    public func minimize(x0: [Double]) -> [Double] {
        return CoordinateDescent(x0: x0)
//        let x0 = Array(repeating: 0.0, count: Settings.N)
//        var result = x0
//        let iters = Hooke(x0: x0, result: &result, rho: 0.1, iterMax: 1000)
//        print("Всего итераций: \(iters)")
//        return result
    }
    
    private func ParabolicMethod(x: [Double], i: Int) -> [Double] {
        let accuracy = 0.01
        var h = 2 * 0.1
        
        var nextX = 0.0
        var currentX = 0.0
        
        var x1 = x
        var x2 = x
        var x3 = x
        
        repeat {
            h /= 2
            
            currentX = x2[i]
            
            x1[i] = currentX - h
            x3[i] = currentX + h
            
            let f1 = CalculateIntegral(w: x1)
            let f2 = CalculateIntegral(w: x2)
            let f3 = CalculateIntegral(w: x3)
            
            let b = (f3 - f1) / (2 * h)
            let a = (f3 - 2 * f2 + f1) / (2 * pow(h, 2))
            
            if (a > 0.001) {
                nextX = currentX - b / (2 * a)
            }
            else {
                let minF = min(f1, f2, f3)
                switch (minF) {
                case f1:
                    nextX = currentX - h
                    break
                case f2:
                    nextX = currentX
                    break
                case f3:
                    nextX = currentX + h
                    break
                default:
                    nextX = 100
                    break
                }
            }
            
            x1[i] = nextX - h
            x2[i] = nextX
            x3[i] = nextX + h
        }
            while (abs(nextX - currentX) > accuracy)
        
        return x2
    }
    
    private func MethodOfTheGoldenRatio(x: [Double], i: Int, a: Double, b: Double, accuracy: Double) -> [Double] {
        
        let phi = (1 + sqrt(5))/2
        
        var _a = a
        var _b = b
        var _x = x
        var x1: Double
        var x2: Double
        var y1, y2: Double
        
        repeat {
            x1 = _b - (_b - _a) / phi
            x2 = _a + (_b - _a) / phi
            
            _x[i] = x1
            y1 = CalculateIntegral(w: _x)
            
            _x[i] = x2
            y2 = CalculateIntegral(w: _x)
            
            if (y1 >= y2) {
                _a = x1
            }
            else {
                _b = x2
            }
        }
            while (abs(_b - _a) > accuracy)
        
        _x[i] = (_a + _b) / 2
        
        return _x
    }
    
    private func CoordinateDescent(x0: [Double]) -> [Double]
    {
        var x = x0
        var counter = 0
        var integral = 0.0
        
        repeat {
            integral = fabs(CalculateIntegral(w: x))
            
            print("Итерация: \(counter)")
            print("Управление: \(x)")
            print("Интеграл: \(integral)")
            
            for i in 0...Settings.N-1 {
                x = ParabolicMethod(x: x, i: i)
//                for i in 0...Settings.N-1 {
//                    if (x[i] >= 10) {
//                        x[i] = 10
//                    }
//                    else if (x[i] <= -10) {
//                        x[i] = -10
//                    }
//                }
                print(i)
            }
            counter += 1
        }
            while (integral > Settings.Accuracy)
        
        return x
    }
    
    private func BestNearBy(delta: [Double], x: inout [Double], prevBest: Double) -> Double {
        var z = x
        var minf = prevBest
        var ftmp = 0.0
        
        for i in 0...Settings.N-1 {
            z[i] = x[i] + delta[i]
            ftmp = CalculateIntegral(w: z)
            if (ftmp < minf) {
                minf = ftmp
            }
            else {
                z[i] = x[i] - delta[i]
                ftmp = CalculateIntegral(w: z)
                if (ftmp < minf) {
                    minf = ftmp
                }
                else {
                    z[i] = x[i]
                }
            }
        }
        
        x = z
        
        print("newf: \(minf)")
        return minf
    }
    
    private func Hooke(x0: [Double], result: inout [Double], rho: Double, iterMax: Int) -> Int {
        var iters = 0
        var iadj = 0
        var stepLength = rho
        var newx = x0
        var xbefore = x0
        var keep = 0
        
        var delta = Array(repeating: rho, count: Settings.N)
        for i in 0...Settings.N-1 {
            if (x0[i] != 0) {
                delta[i] *= x0[i]
            }
        }
        
        var fbefore = CalculateIntegral(w: newx)
        var newf = fbefore
        
        print("Итерация: \(iters)")
        print("Управление: \(newx)")
        print("Интеграл: \(fbefore)")
        
        while (newf > Settings.Accuracy) {

            iters += 1
            iadj += 1
            
            newx = xbefore
            newf = BestNearBy(delta: delta, x: &newx, prevBest: fbefore)

            keep = 1
            while (newf < fbefore && keep == 1) {
                iadj = 0
                for i in 0...Settings.N-1 {
                    if (newx[i] <= xbefore[i]) {
                        delta[i] = -abs(delta[i])
                    }
                    else {
                        delta[i] = abs(delta[i])
                    }
                    
                    let tmp = xbefore[i]
                    xbefore[i] = newx[i]
                    newx[i] = newx[i] + newx[i] - tmp
                }
                fbefore = newf
                newf = BestNearBy(delta: delta, x: &newx, prevBest: fbefore)
                
                if (newf >= fbefore) {
                    break
                }
                
                keep = 0
                for i in 0...Settings.N-1 {
                    keep = 1
                    
                    let absv = abs(newx[i] - xbefore[i])
                    if (absv > 0.5 * abs(delta[i])) {
                        break
                    }
                    else {
                        keep = 0
                    }
                }
            }
            
            if (stepLength >= Settings.Accuracy / 10000 && newf >= fbefore) {
                stepLength *= rho
                for i in 0...Settings.N-1 {
                    delta[i] *= rho
                }
            }
            
            print("Итерация: \(iters)")
            print("Управление: \(xbefore)")
            print("Интеграл: \(fbefore)")
        }
        
        result = xbefore
        return iters
    }
    
    public func CalculateIntegral(w: [Double]) -> Double {
        var result = 0.0
        
        let solution = Solve(w: w)
        
        for k in 1...K {
            for m in 0...M-1 {
                let a1 = pow(solution[N-1][k][m], 2)
                let aa1 = (solution[N-1][k][m] - solution[N-2][k][m])
                let a2 = pow(aa1 / tau, 2)

                let b = (Double(m) + 0.5) * pow(h_r, 2) * h_phi

                result += (a1 + a2) * b
            }
        }
        
        return result
    }
}
