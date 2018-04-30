import Foundation

public class Worker {
    public init() {
        let w0 = Array(repeating: 0.0, count: Settings.N)
        let equation = Equation(_T: Settings.T, _K: Settings.K, _M: Settings.M, _N: Settings.N, _R: Settings.R, _h_r: Settings.h_r(), _h_phi: Settings.h_phi(), _tau: Settings.tau())
        let solution = equation.Solve(w: w0)
        
        var minimizedSolution = solution
        var w = w0
        var workTimeSec = 0.0
        
        if (Settings.Minimize) {
            let start = Date()
            w = equation.minimize(x0: w0)
            minimizedSolution = equation.Solve(w: w)
            workTimeSec = Date().timeIntervalSince(start)
        }
        
        let shared = Shared(
            _T: Settings.T,
            _K: Settings.K,
            _M: Settings.M,
            _N: Settings.N,
            _R: Settings.R,
            _h_r: Settings.h_r(),
            _h_phi: Settings.h_phi(),
            _tau: Settings.tau(),
            _w0: w0,
            _w: w,
            _soultion: solution,
            _minimizedSolution: minimizedSolution,
            _workTimeSec: workTimeSec,
            _actuatorIndex: Settings.ActuatorIndex,
            _accuracy: Settings.Accuracy,
            _exampleNumber: Settings.ExampleNumber,
            _actuatorType: Settings.ActuatorType);
        
        let payload =  String(data: try! JSONEncoder().encode(shared), encoding: String.Encoding.utf8) as String!
        
        let formatter = DateFormatter()
        formatter.dateFormat = "HH-mm-ss dd.MM.yyyy"
        let fileName = formatter.string(from: Date())
        FileHelper.AppendFile(fileName: fileName + ".txt", text: payload!)
    }
}
