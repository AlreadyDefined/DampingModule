//
//  Settings.swift
//  DampingModule
//
//  Created by Irina Filippova.
//  Copyright © 2017 Irina Filippova. All rights reserved.
//

import Foundation

public class Settings {
    public static var T = 2.0
    public static var K = 3
    public static var M = 5
    public static var N = 800
    
    public static var R = 1.0
    
    public static var ExampleNumber = 7
    
    public static var Minimize = true
    
    public static var ActuatorIndex = 1
    
    public static var Accuracy = 0.01
    
    public static var ActuatorType = ActType.point
    
    public static func h_r() -> Double {
        return Double(R) / (Double(M) - 0.5)
    }
    
    public static func h_phi() -> Double {
        return .pi / Double(K-1)
    }
    
    public static func tau() -> Double {
        return T / Double(N-1)
    }
    
    private static func Example1(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        let t = n * Settings.tau()
        
        switch(type) {
        case .h0:
            return 0.0
        case .h1:
            return r * (1 - r)
        case .f:
            return 2*t - (1 - 2 * r) * t / r
        case .u:
            return r * (1 - r) * t
        }
    }
    
    private static func Example2(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        let t = n * Settings.tau()
        
        switch(type) {
        case .h0:
            return sin(.pi * r)
        case .h1:
            return sin(.pi * r)
        case .f:
            return (pow(.pi, 2) * sin(.pi * r) - .pi /
                r * cos(.pi * r)) * (t + 1)
        case .u:
            return (t + 1) * sin(.pi * r)
        }
    }
    
    private static func Example3(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        let t = n * Settings.tau()
        
        switch(type) {
        case .h0:
            return 0.0
        case .h1:
            return 1 - pow(r, 2)
        case .f:
            return 4 * t
        case .u:
            return t * (1 - pow(r, 2))
        }
    }
    
    private static func Example4(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        let t = n * Settings.tau()
        
        switch(type) {
        case .h0:
            return 1.0/8.0 * cos(.pi * r / 2.0)
        case .h1:
            return 0
        case .f:
            let const = cos(t) * cos(.pi * r / 2.0)
            return -1.0/8.0 * const +
                pow(.pi, 2) / 32.0 * const +
                .pi/16.0 * cos(t) * sin(.pi * r / 2.0) / r
        case .u:
            return 1.0/8.0 * cos(t) * cos(.pi * r / 2.0)
        }
    }
    
    private static func Example5(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        
        switch(type) {
        case .h0:
            return r <= 1.0/3.0 ? 1 : 3.0/2.0 * (1 - r)
        case .h1:
            return 1 - r
        default:
            return 0
        }
    }
    
    private static func Example6(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        
        switch(type) {
        case .h0:
            return 1.0 / 2.0 * sin(.pi * r)
        case .h1:
            return 0
        default:
            return 0
        }
    }
    
    private static func Example7(type: FunctionType,
                                 m: Int, n: Double) -> Double {
        let r = (Double(m) + 0.5) * Settings.h_r()
        
        switch(type) {
        case .h0:
            return 1.0 / 2.0 * (1 - r * r)
        case .h1:
            return 0
        default:
            return 0
        }
    }
    
    public static func Example(type: FunctionType,
                               m: Int, n: Double = 0.0) -> Double {
        switch(Settings.ExampleNumber) {
        case 1:
            return Settings.Example1(type: type, m: m, n: n)
        case 2:
            return Settings.Example2(type: type, m: m, n: n)
        case 3:
            return Settings.Example3(type: type, m: m, n: n)
        case 4:
            return Settings.Example4(type: type, m: m, n: n)
        case 5:
            return Settings.Example5(type: type, m: m, n: n)
        case 6:
            return Settings.Example6(type: type, m: m, n: n)
        case 7:
            return Settings.Example7(type: type, m: m, n: n)
        default:
            return 0.0
        }
    }
    
    public enum FunctionType {
        case h0
        case h1
        case f
        case u
    }
    
    public enum ActType {
        case none
        case circular
        case point
    }
}
