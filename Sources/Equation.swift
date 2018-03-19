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
    
    private func calculate_f(w: [Double]) -> Array<Array<[Double]>> {
        var f = Array(repeating: Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3), count: N + 1)
        
        for n in 0...N {
            for k in 1...K+1 {
                //for m in 0...M {
                f[n][k][Settings.ActuatorIndex] = w[n]
                //Вариант для сравнения с аналитическим значением
                //f[n][k][m] = Settings.Example(type: Settings.FunctionType.f, m: m, n: Double(n))
                //}
            }
        }
        
        return f
    }
    
    private func calculate_h0() -> Array<[Double]> {
        var h0 = Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3)
        
        for k in 1...K+1 {
            for m in 0...M {
                h0[k][m] = Settings.Example(type: Settings.FunctionType.h0, m: m)
            }
        }
        
        return h0
    }
    
    private func calculate_h1() -> Array<[Double]> {
        var h1 = Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3)
        
        for k in 1...K+1 {
            for m in 0...M {
                h1[k][m] = Settings.Example(type: Settings.FunctionType.h1, m: m)
            }
        }
        
        return h1
    }
    
    private func getFirstTimeLayerOld(prev: Array<[Double]>, f: Array<[Double]>, h1: Array<[Double]>) -> Array<[Double]> {
        var next = Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3)
        // next[k][M] остаётся равным нулю
        
        for k in 1...K+1 {
            //u_k0^1
            next[k][0] =    prev[k][0] +
                tau * h1[k][0] +
                pow(tau / h_r, 2) * (prev[k][1] - prev[k][0]) +
                2 * pow(tau / (h_r * h_phi), 2) * (prev[k+1][0] - 2 * prev[k][0] + prev[k-1][0]) + pow(tau, 2) * f[k][0] * 0.5
            
            //u_km^1
            for m in 1...M-1 {
                next[k][m] =    prev[k][m] +
                    tau * h1[k][m] +
                    0.5 * pow(tau / h_r, 2) * bracesOld(prev: prev, f: f, k: k, m: m)
            }
        }
        
        return next
    }
    
    private func getFirstTimeLayerNew(prev: Array<[Double]>, f: Array<[Double]>, h1: Array<[Double]>) -> Array<[Double]> {
        var next = Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3)
        // next[k][M] остаётся равным нулю
        
        for k in 1...K+1 {
            //u_km^1
            for m in 0...M-1 {
                next[k][m] =    prev[k][m] +
                    tau * h1[k][m] +
                    0.5 * pow(tau / h_r, 2) * bracesNew(prev: prev, f: f, k: k, m: m)
            }
        }
        
        return next
    }
    
    private func getOtherTimeLayersOld(cur: Array<[Double]>, prev: Array<[Double]>, f: Array<[Double]>) -> Array<[Double]> {
        // next[k][M] остаётся равным нулю
        var next = Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3)
        
        //u_k0^n+1
        for k in 1...K+1 {
            next[k][0] =    2 * cur[k][0] -
                prev[k][0] +
                2 * pow(tau / h_r, 2) * (cur[k][1] - cur[k][0]) +
                4 * pow(tau / (h_r * h_phi), 2) * (cur[k+1][0] - 2 * cur[k][0] + cur[k-1][0]) +
                pow(tau, 2) * f[k][0]
        }
        
        //u_km^n+1
        for k in 1...K+1 {
            for m in 1...M-1 {
                next[k][m] =    2 * cur[k][m] -
                    prev[k][m] +
                    pow(tau / h_r, 2) * bracesOld(prev: prev, f: f, k: k, m: m)
            }
        }
        
        return next
    }
    
    private func getOtherTimeLayersNew(cur: Array<[Double]>, prev: Array<[Double]>, f: Array<[Double]>) -> Array<[Double]> {
        // next[k][M] остаётся равным нулю
        var next = Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3)
        
        //u_km^n+1
        for k in 1...K+1 {
            for m in 0...M-1 {
                next[k][m] =    2 * cur[k][m] -
                    prev[k][m] +
                    pow(tau / h_r, 2) * bracesNew(prev: prev, f: f, k: k, m: m)
            }
        }
        
        return next
    }
    
    private func bracesOld(prev: Array<[Double]>, f: Array<[Double]>, k: Int, m: Int) -> Double {
        
        let const1 = prev[k][m+1] - 2 * prev[k][m] + prev[k][m - 1]
        let const2 = (prev[k][m + 1] - prev[k][m - 1]) / (2 * Double(m) + 1)
        let const3 = (prev[k + 1][m] - 2 * prev[k][m] + prev[k - 1][m]) / pow((Double(m) + 0.5) * h_phi, 2)
        
        return  (const1 +
            const2 +
            const3) +
            f[k][m] * pow(h_r, 2)
    }
    
    private func bracesNew(prev: Array<[Double]>, f: Array<[Double]>, k: Int, m: Int) -> Double {
        
        let temp1 = (Double(m) + 1) * (prev[k][m+1] - prev[k][m])
        let temp2 = m == 0 ? 0 : Double(m) * (prev[k][m] - prev[k][m - 1])
        let const1 = (temp1 - temp2) / (Double(m) + 0.5)
        let const2 = (prev[k + 1][m] - 2 * prev[k][m] + prev[k - 1][m]) / pow((Double(m) + 0.5) * h_phi, 2)
        
        return  (const1 +
            const2) +
            f[k][m] * pow(h_r, 2)
    }
    
    public func solveOld(w: [Double]) -> Array<Array<[Double]>> {
        var solution = Array(repeating: Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3), count: N+1)
        
        let f = calculate_f(w: w)
        let h1 = calculate_h1()
        
        for n in 0...N {
            switch(n) {
            case 0:
                solution[n] = calculate_h0()
                break;
            case 1:
                solution[n] = getFirstTimeLayerOld(prev: solution[0], f: f[0], h1: h1)
                break;
            default:
                solution[n] = getOtherTimeLayersOld(cur: solution[n - 1], prev: solution[n - 2], f: f[n - 1])
                break;
            }
            
            shift(u: &solution[n])
        }
        
        
        //print("Решение: \(solution)")
        //let maxDiff = MathHelper.calculateMaxDiff(actual: solution, expected: MathHelper.calculateExactFunction())
        //print("exact: \(MathHelper.calculateExactFunction()[N])")
        //print("diff: \(maxDiff)")
        return solution
    }
    
    public func solveNew(w: [Double]) -> Array<Array<[Double]>> {
        var solution = Array(repeating: Array(repeating: Array(repeating: 0.0, count: M+1), count: K+3), count: N+1)
        
        let f = calculate_f(w: w)
        let h1 = calculate_h1()
        
        for n in 0...N {
            switch(n) {
            case 0:
                solution[n] = calculate_h0()
                break;
            case 1:
                solution[n] = getFirstTimeLayerNew(prev: solution[0], f: f[0], h1: h1)
                break;
            default:
                solution[n] = getOtherTimeLayersNew(cur: solution[n - 1], prev: solution[n - 2], f: f[n - 1])
                break;
            }
            
            shift(u: &solution[n])
        }
        
        
        //print("Решение: \(solution)")
        //let maxDiff = MathHelper.calculateMaxDiff(actual: solution, expected: MathHelper.calculateExactFunction())
        //print("exact: \(MathHelper.calculateExactFunction()[N])")
        //print("diff: \(maxDiff)")
        return solution
    }
    
    private func shift(u: inout Array<[Double]>) {
        for m in 0...M {
            u[0][m] = u[2][m]
            u[K + 2][m] = u[K][m]
        }
    }
    
    public func minimize(x0: [Double]) -> [Double] {
        return CoordinateDescent(N: N, x0: x0)
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
            
            let f1 = calculateIntegral(w: x1)
            let f2 = calculateIntegral(w: x2)
            let f3 = calculateIntegral(w: x3)
            
            let b = (f3 - f1) / (2 * h)
            let a = (f3 - 2 * f2 + f1) / (2 * pow(h, 2))
            
            if (a > 0.01) {
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
            y1 = calculateIntegral(w: _x)
            
            _x[i] = x2
            y2 = calculateIntegral(w: _x)
            
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
    
    private func CoordinateDescent(N: Int, x0: [Double]) -> [Double]
        //minimizes N-dimensional function f; x0 - start point
    {
        var x = x0
        var counter = 0
        var integral = 0.0
        
        repeat {
            integral = fabs(calculateIntegral(w: x))
            
            print("Итерация: \(counter)")
            print("Управление: \(x)")
            print("Интеграл: \(integral)")
            
            for i in 0...N {
                //x = MethodOfTheGoldenRatio(x: x, i: i, a: -10, b: 10, accuracy: accuracy)
                x = ParabolicMethod(x: x, i: i)
                print(i)
            }
            counter += 1
        }
            while (integral > Settings.Accuracy)
        
        return x
    }
    
    public func calculateIntegral(w: [Double]) -> Double {
        var result = 0.0
        
        let solution = //Settings.NewAlgorithm
            //? solveNew(w: w) :
            solveNew(w: w)
        
        for k in 1...K+1 {
            for m in 1...M { // m in 1...M-1
                let a1 = pow(solution[N][k][m], 2)
                let aa1 = (solution[N][k][m] - solution[N-1][k][m])
                let a2 = pow(aa1 / tau, 2)
                
                let b = (Double(m) + 0.5) * h_r * h_phi
                
                result += (a1 + a2) * b
                
                //                    let a1 = (solution[N][k][m+1] - solution[N][k][m-1]) / (2 * h_r)
                //                    let a2 = (solution[N][k][m] - solution[N-1][k][m]) / tau
                //
                //                    result += (pow(a1, 2) + pow(a2, 2))
            }
        }
        
        return result
    }
}
