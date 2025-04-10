using MeshGrid








function main()
    # Test the meshgrid function
    x = [1, 2, 3]
    y = [4, 5, 6]
    X, Y = meshgrid(x, y)
    println(X)
    println(Y)
end





if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
    

