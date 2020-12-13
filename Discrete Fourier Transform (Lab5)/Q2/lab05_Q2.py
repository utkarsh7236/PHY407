import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft2, irfft2

############## Q7.9 (a) ##############
# read the txt file into 2D array
image = np.loadtxt("blur.txt")

# plot the image
plt.figure(figsize=(6, 6))
plt.imshow(image, cmap='gray')
plt.title('Blurred image', fontsize=16)
plt.xlabel('X pixel index', fontsize=14)
plt.ylabel('Y pixel index', fontsize=14)
plt.tight_layout()
plt.savefig('Q2a.pdf')


############## Q7.9 (b) ##############
# get image size
height, width = image.shape

# calculate point spread function
# set sigma
sigma = 25
# initialize empty array for psf
psf = np.empty([height, width]) 
# loop through width and height and calculate psf. Code from the handout
for i in range(height):
    ival = i
    if ival > height // 2:
        ival = height - ival # flip the sign for the bottom half
    for j in range(width):
        jval = j
        if jval > width//2:
            jval = width - jval # flip the sign of the right half
        # calculate and store gaussian
        psf[i,j] = np.exp(- (ival**2 + jval**2) / (2 * sigma**2))

# visualize gaussian kernel
plt.figure(figsize=(6, 6))
plt.imshow(psf, cmap='gray', vmax=0.2)
plt.title('Density plot (Gussian with $\sigma=25$)', fontsize=16)
plt.xlabel('X pixel index', fontsize=14)
plt.ylabel('Y pixel index', fontsize=14)
plt.tight_layout()
plt.colorbar()
plt.savefig('Q2b.pdf')


############## Q7.9 (c) ##############
# Fourier transform image and psf
ft_image = rfft2(image)
ft_psf = rfft2(psf)
# change values less that 1e-3 in ft_psf to 1 to avoid division by 0
thresh = 1e-3
ft_psf[ft_psf <= thresh] = 1
# divide ft_image by ft_psf
div_result = ft_image / (ft_psf * height * width)

# perform inverse fourier transform to get unblurred image
unblurred_image = irfft2(div_result)

# plot unblurred image
plt.figure(figsize=(6, 6))
plt.imshow(unblurred_image, cmap='gray')
plt.title('Unblurred image', fontsize=16)
plt.xlabel('X pixel index', fontsize=14)
plt.ylabel('Y pixel index', fontsize=14)
plt.tight_layout()
plt.savefig('Q2c.pdf')

