// swift-tools-version:4.0

import PackageDescription

let package = Package(
    name: "DampingModule",
    products: [
        .library(name: "DampingModule", targets: ["DampingModule"]),
        ],
    targets: [
        .target(
            name: "DampingModule",
            dependencies: [],
            path: "Sources"
        ),
    ]
)


