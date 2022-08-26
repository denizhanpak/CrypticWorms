membership = [1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2]
membership = 1:9
nodecolor = [colorant"lightseagreen", colorant"orange"]
colors = [colorant"blue",colorant"red",colorant"green",
    colorant"yellow",colorant"purple",colorant"pink",
    colorant"orange",colorant"brown",colorant"white"]
# membership color
nodefillc = colors[membership]
gplot(g, nodefillc=nodefillc,layout=grid_layout)
