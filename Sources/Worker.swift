import Foundation

public class Worker {
    public init() {
        let w0 = Array(repeating: 0.0, count: Settings.N + 1)
        var w = Array(repeating: 0.0, count: Settings.N + 1)
        
        let equation = Equation(_T: Settings.T, _K: Settings.K, _M: Settings.M, _N: Settings.N, _R: Settings.R, _h_r: Settings.h_r(), _h_phi: Settings.h_phi(), _tau: Settings.tau())
        
        let solution = Settings.NewAlgorithm
            ? equation.solveNew(w: w0)
            : equation.solveOld(w: w0)
        
        let minimized = equation.minimize(x0: w0)
        //let minimizedSolution = solution
        
        w = minimized
        let minimizedSolution = Settings.NewAlgorithm
            ? equation.solveNew(w: w)
            : equation.solveOld(w: w)
        
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
            _newAlgorithm: Settings.NewAlgorithm,
            _actuatorIndex: Settings.ActuatorIndex,
            _accuracy: Settings.Accuracy);
        
        let payload =  String(data: try! JSONEncoder().encode(shared), encoding: String.Encoding.utf8) as String!
        
        let formatter = DateFormatter()
        formatter.dateFormat = "HH/mm dd.MM.yyyy"
        let fileName = formatter.string(from: Date())
        FileHelper.AppendFile(fileName: fileName + ".txt", text: payload!)
    }
}
