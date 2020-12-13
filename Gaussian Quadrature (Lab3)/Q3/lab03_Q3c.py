#Author: Aslesha, Utkarsh

import numpy as np
import struct
import matplotlib.pyplot as plt


# defining our image function for intensity
def I(theta, dwdx, dwdy):
    """
    Return intensity using the equation in pg. 212 of the textbook.
    """
    sq = np.sqrt((dwdx)**2 + (dwdy)**2 + 1)
    return - ((np.cos(theta) * dwdx + np.sin(theta) * dwdy) / (sq))


# Defining our initial read function
def main():
    # establish file path
    file_name = "N46E006.hgt"

    # create an empty array for w
    w = []

    # open file and initialise for loop
    with open(file_name, "rb") as f:

        # loop over each column
        for i in range(1201):
            row = []

            # loop over each row
            for j in range(1201):

                # read byte at specific pixel location 
                buf = f.read(2)

                # convert bytes to binary data and append 
                row.append(struct.unpack('>h', buf)[0])
            w.append(row)

    # show the tile
    plt.figure(figsize=(10,8))
    plt.imshow(w, origin = 'upper',  extent=[6, 7, 46, 47], vmin = 0.1)
    plt.ylabel("Latitude $^{\circ}N$")
    plt.xlabel("Longitude $^{\circ}E$")
    plt.title("Elevation tile N46E006")
    plt.colorbar()
    plt.savefig("Q3a.pdf")
    

    #empty array for dwdx
    dwdx = []

    # h step defined 
    h = 420

    # calculate the gradient with respect to dx
    for i in range(len(w)):

        # creating empty array for the row gradients 
        row_grad = []

        # for loop over the rows
        for j in range(len(w[0])):

            # using forward difference if j is the first item
            if j == 0: 
                row_grad.append((w[i][j+1] - w[i][j]) / h)
                
            # use backwards difference is j is the last item
            elif j == len(w[0]) - 1: 
                row_grad.append((w[i][j] - w[i][j-1]) / h)
            
            # central difference otherwise
            else:
                row_grad.append((w[i][j+1] - w[i][j-1]) / 2 * h)

        # append our row onto the partial derivative
        dwdx.append(row_grad)

    # repeat steps with dwdy
    dwdy = []
    
    # define our h 
    h = 420

    # calculate the gradient with respect to dy
    for i in range(len(w)):

        # creating empty array for the column gradients 
        col_grad = []

        # for loop over the columns        
        for j in range(len(w[0])):

            # using forward difference if i is the first item
            if i == 0:
                col_grad.append((w[i+1][1] - w[i][j]) / h)

            # use backwards difference is i is the last item
            elif i == len(w[0]) - 1:
                col_grad.append((w[i][j] - w[i-1][j]) / h)

            # central difference otherwise    
            else:
                col_grad.append((w[i+1][j] - w[i-1][j]) / 2 * h)

        # append our column onto the partial derivative
        dwdy.append(col_grad)

    # set value given for theta 
    theta = - 5 * np.pi / 6
    
    # use intensity function to calculate image for intensity
    intensity = I(theta, np.array(dwdx), np.array(dwdy))

    # plot intensity
    plt.figure(figsize=(10,8))
    plt.imshow(intensity, cmap='gray', origin = 'upper',  extent=[6, 7, 46, 47])
    plt.ylabel("Latitude $^{\circ}N$")
    plt.xlabel("Longitude $^{\circ}E$")
    plt.title("Intensity map")
    plt.colorbar()
    plt.savefig("Q3b.pdf")
    

if __name__ == "__main__":
    main() # running our experiment