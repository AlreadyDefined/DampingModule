//
//  Shared.swift
//  HeatConductionEquation
//
//  Created by Ирина Филиппова on 13.02.2018.
//  Copyright © 2018 Ирина Филиппова. All rights reserved.
//

import Foundation

public class Shared : Codable {
    init(_T: Double = 0.0,
         _K: Int = 0,
         _M: Int = 0,
         _N: Int = 0,
         _R: Double = 0,
         _h_r: Double = 0.0,
         _h_phi: Double = 0.0,
         _tau: Double = 0.0,
         _w0: Array<Double> = [],
         _w: Array<Double> = [],
         _soultion: Array<Array<[Double]>> = [[[]]],
         _minimizedSolution: Array<Array<[Double]>> = [[[]]],
         _workTimeSec: Double = 0.0,
         _actuatorIndex: Int = 0,
         _accuracy: Double = 0.0,
         _exampleNumber: Int = 0,
         _actuatorType: Settings.ActType = Settings.ActType.none,
         _iters: Int = 0) {
        T = _T
        K = _K
        M = _M
        N = _N
        
        R = _R
        
        h_r = _h_r
        h_phi = _h_phi
        tau = _tau
        
        w0 = _w0
        w = _w
        
        solution = _soultion
        minimizedSolution = _minimizedSolution
        
        workTimeSec = _workTimeSec
        actuatorIndex = _actuatorIndex
        accuracy = _accuracy
        exampleNumber = _exampleNumber
        actuatorType = String(describing: _actuatorType)
        iters = _iters
    }
    
    public static var instance = Shared()
    
    public var T: Double
    public var K: Int
    public var M: Int
    public var N: Int
    
    public var R: Double
    
    public var h_r: Double
    public var h_phi: Double
    public var tau: Double
    
    public var w0 = Array<Double> ()
    public var w = Array<Double> ()
    
    public var solution = Array<Array<[Double]>>()
    public var minimizedSolution = Array<Array<[Double]>>()

    public var workTimeSec: Double
    public var actuatorIndex: Int
    public var accuracy: Double
    public var exampleNumber: Int
    public var actuatorType: String
    public var iters: Int
}

