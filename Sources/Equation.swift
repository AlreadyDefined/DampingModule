import Foundation

class Equation {
    var T: Double = 0.0
    var K: Int = 0
    var M: Int = 0
    var N: Int = 0
    var R: Double = 0.0
    
    var h_r: Double = 0.0
    var h_phi: Double = 0.0
    var tau: Double = 0.0
    var h0: Array<[Double]> = []
    var h1: Array<[Double]> = []
    var iters: Int = 0
    
    init(_T: Double, _K: Int, _M: Int, _N: Int, _R: Double, _h_r: Double, _h_phi: Double, _tau: Double){
        T = _T
        K = _K
        M = _M
        N = _N
        R = _R
        
        h_r = _h_r
        h_phi = _h_phi
        tau = _tau
        
        (h0, h1) = self.initialize()
    }
    
    private func calculate_f(w: [Double]) -> Array<Array<[Double]>> {
        var f = Array(repeating: Array(repeating: Array(repeating: 0.0, count: M), count: K+2), count: N-1)
        
//        switch (Settings.ActuatorType) {
//        case Settings.ActType.none:
//            for n in 0...N-1 {
//                for k in 1...K {
//                    for m in 0...M-1 {
//                        f[n][k][m] = Settings.Example(type: Settings.FunctionType.f, m: m, n: Double(n))
//                    }
//                }
//            }
//            break
//        case Settings.ActType.circular:
//            for n in 0...N-1 {
//                for k in 1...K {
//                    f[n][k][Settings.ActuatorIndex] = w[n]
//                }
//            }
//            break
//        case Settings.ActType.point:
            for n in 0...N-2 {
//                f[n][2][Settings.ActuatorIndex] = w[n]
                f[n][1][Settings.ActuatorIndex] = w[n]
//                f[n][0][Settings.ActuatorIndex] = w[n]
//
//                f[n][2][Settings.ActuatorIndex + 1] = w[n]
//                f[n][1][Settings.ActuatorIndex + 1] = w[n]
//                f[n][0][Settings.ActuatorIndex + 1] = w[n]
//
//                f[n][2][Settings.ActuatorIndex + 2] = w[n]
//                f[n][1][Settings.ActuatorIndex + 2] = w[n]
//                f[n][0][Settings.ActuatorIndex + 2] = w[n]
            }
//        }

        return f
    }
    
    private func initialize() -> (h0: Array<[Double]>, h1: Array<[Double]>){
        var h0 = Array(repeating: Array(repeating: 0.0, count: M), count: K+2)
        var h1 = h0

        for k in 1...K {
            for m in 0...M-1 {
                h0[k][m] = Settings.Example(type: Settings.FunctionType.h0, m: m)
                h1[k][m] = Settings.Example(type: Settings.FunctionType.h1, m: m)
            }
        }
        
        return (h0, h1)
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
        let f = calculate_f(w: w)
        
        for n in 0...N-1 {
            switch(n) {
            case 0:
                solution[n] = h0
                break
            case 1:
                solution[n] = GetFirstTimeLayer(prev: solution[0], f: f[0], h1: h1)
                break
            default:
                solution[n] = GetOtherTimeLayers(cur: solution[n - 1], prev: solution[n - 2], f: f[n - 1])
                break
            }
            
            Shift(u: &solution[n])
        }
        
//        let maxDiff = MathHelper.calculateMaxDiff(actual: solution[N-1], expected: MathHelper.calculateExactFunction())
//        print("actual: \(solution[N-1])")
//        print("exact: \(MathHelper.calculateExactFunction())")
//        print("diff: \(maxDiff)")
//        let maxDiff = MathHelper.calculateMaxDiff1(actual: solution[N-1], expected: MathHelper.calculateExactFunction())
//        print("exact: \(MathHelper.calculateExactFunction())")
//        print("diff: \(maxDiff)")
        
//        for n in 0...N-1 {
//            for m in 0...M-1 {
//                for k in 0...K+1 {
//                    print("u_\(k) = \(solution[n][k][m])")
//                }
//                print("m = \(m)")
//            }
//            print("n = \(n)")
//
//        }
        
        return solution
    }
    
    private func Shift(u: inout Array<[Double]>) {
        for m in 0...M-1 {
            u[0][m] = u[2][m]
            u[K + 1][m] = u[K - 1][m]
        }
    }
    
    public func Minimize(x0: [Double]) -> [Double] {
        return CoordinateDescent(x0: x0)
    }
    
    private func ParabolicMethod(x: [Double], i: Int) -> Double {
        let accuracy = 0.01
        var h = 2.0

        var currentX = 0.0
        
        var x1 = x
        var x2 = x
        var x3 = x
        
        repeat {
            currentX = x2[i]
            
            x1[i] = currentX - h
            x3[i] = currentX + h

            let f1 = CalculateIntegral(w: x1)
            let f2 = CalculateIntegral(w: x2)
            let f3 = CalculateIntegral(w: x3)

            let b = (f3 - f1) / (2 * h)
            let a = (f3 - 2 * f2 + f1) / (2 * pow(h, 2))
            
            if (a > 0.001) {
                x2[i] = currentX - b / (2 * a)
            }
            else {
                let minF = min(f1, f2, f3)
                if (minF == f1) {
                    x2[i] = currentX - h
                }
                else if (minF == f2) {
                    x2[i] = currentX
                }
                else {
                    x2[i] = currentX + h
                }
            }
            
            if (abs(x2[i] - currentX) < accuracy) {
                h /= 2
            }
            
            if (h < accuracy) {
                break
            }
        }
        while (true)

        return x2[i]
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
            
            for i in 0...N-2 {
                x[i] = ParabolicMethod(x: x, i: i)
                print(i)
            }
            counter += 1
        }
        while (integral > Settings.Accuracy)
        iters = counter
        return x
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
