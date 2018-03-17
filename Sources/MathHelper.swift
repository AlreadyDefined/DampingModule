//
//  Math.swift
//  HeatConductionEquation
//
//  Created by Ирина Филиппова on 04.10.17.
//  Copyright © 2017 Ирина Филиппова. All rights reserved.
//

import Foundation
class MathHelper {
    public static func calculateMaxDiff(actual: Array<Array<[Double]>>, expected: Array<Array<[Double]>>) -> Double {
        var diff = Array(repeating: Array(repeating: Array(repeating: 0.0, count: Settings.M), count: Settings.K + 2), count: Settings.N + 1)
        
        var max = 0.0
        
        for n in 0...Settings.N {
            for k in 1...Settings.K {
                for m in 0...Settings.M-1 {
                    diff[n][k][m] = abs(expected[n][k][m] - actual[n][k][m])
                    if (diff[n][k][m] > max) {
                        max = diff[n][k][m]
                    }
                }
            }
        }
        
        return max
    }

    public static func calculateExactFunction() -> Array<Array<[Double]>> {
        var result = Array(repeating: Array(repeating: Array(repeating: 0.0, count: Settings.M), count: Settings.K + 2), count: Settings.N + 1)
        
        for n in 0...Settings.N {
            for k in 0...Settings.K {
                for m in 0...Settings.M-1 {
                    result[n][k][m] = Settings.Example(type: Settings.FunctionType.u, m: m, n: Double(n))
                }
            }
        }
        
        return result
    }
}
