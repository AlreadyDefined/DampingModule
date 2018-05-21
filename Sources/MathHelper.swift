//
//  Math.swift
//  DampingModule
//
//  Created by Irina Filippova.
//  Copyright © 2017 Irina Filippova. All rights reserved.
//

import Foundation
class MathHelper {
    public static func calculateMaxDiff(actual: Array<[Double]>,
        expected: Array<[Double]>) -> Double {
        
        var maxDiff = 0.0
        
        for k in 1...Settings.K {
            for m in 0...Settings.M-1 {
                maxDiff = max(abs(expected[k][m] - actual[k][m]), maxDiff)
            }
        }
        
        return maxDiff
    }
    
    public static func CheckRungeRule(actual: Array<[Double]>,
        expected: Array<[Double]>) -> Double {
        
        let diff1 = expected[1][10]
        let diff2 = actual[1][10]
        print("expected = \(diff1), actual = \(diff2)")
        return abs(diff1 - diff2)
    }

    public static func CalculateExactFunction() -> Array<[Double]> {
        var result = Array(repeating: Array(repeating: 0.0,
            count: Settings.M), count: Settings.K+2)
        
        for k in 1...Settings.K {
            for m in 0...Settings.M-1 {
                result[k][m] = Settings.Example(type: Settings.FunctionType.u,
                     m: m, n: Double(Settings.N - 1))
            }
        }
        
        return result
    }
}
