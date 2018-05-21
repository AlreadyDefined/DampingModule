//
//  FileHelper.swift
//  DampingModule
//
//  Created by Irina Filippova.
//  Copyright © 2017 Irina Filippova. All rights reserved.
//

import Foundation

public class FileHelper {
    
    public static func ReadFile(url: URL) -> String {
        var result = ""

        do {
            result = try String(contentsOf: url, encoding: .utf8)
        }
        catch {
            print("Error handling \(url)")
        }
        
        return result
    }
    
    public static func AppendFile(fileName: String, text: String) {
        
        let fileURL = URL(fileURLWithPath: fileName)
        
        
        if let fileHandle = FileHandle(forWritingAtPath: fileURL.path) {
            fileHandle.write(text.data(using: .utf8)!)
        }
        else {
            do {
                try text.write(to: fileURL, atomically: true, encoding: .utf8)
            } catch {
                print("Error creating \(fileURL)")
            }
        }
    }
}
