import numpy as np
from PIL import Image

def ImageToFFT(image):
    """
    image: a 2D numpy array
    Return: a 2D numpy array, as the FFT of the image

    You can use np.fft API.
    """
    # TODO: your code here
    # raise NotImplementedError

    img = np.fft.fft2(image)
    img = np.fft.fftshift(img)

    spectrum = np.log(1 + np.abs(img))

    return spectrum


def SampleByThreshold(size, threshold=0.8):
    """
    size: shape of the image
    threshold: threshold to set the value
    Return: an image of white noise
    """
    # TODO: your code here
    # raise NotImplementedError

    img = np.random.rand(*size)

    img[img < threshold] = 0
    img[img >= threshold] = 255

    return img.astype(np.uint8)

def SampleByWindow(size, window_size=(4, 4)):
    """
    size: shape of the image
    window_size: shape of subgrid
    Return: an image of white noise
    """
    # TODO: your code here
    # raise NotImplementedError

    img = np.zeros(size, np.uint8)

    H, W = size
    dh, dw = window_size

    for i in range(0, H, dh):
        for j in range(0, W, dw):
            idx = np.random.randint(0, dh*dw)

            img[i + idx//dw, j + idx%dw] = 255

    return img

def PoissonDiskSampling(size, radius, max_iter=100):
    """
    size: shape of the image
    radius: maximum distance between points
    max_iter: maximum number of iterations to try to add one point
    Reutrn: an image of sampled points
    """
    # TODO: your code here
    # raise NotImplementedError

    def norm2(dx, dy):
        """A helper function to find L2 norm^2.
        """

        return dx*dx + dy*dy

    H, W = size

    ## grid initialization
    cell_size = radius / np.sqrt(2)
    GX = int(np.ceil(H/cell_size))
    GY = int(np.ceil(W/cell_size))
    grid = np.zeros((GX, GY), dtype=np.int_) - 1

    ## select the first sample
    x = np.random.randint(0, H)
    y = np.random.randint(0, W)

    count = 0

    grid[int(x / cell_size), int(y / cell_size)] = count

    active_list = [(x, y)]
    sample_list = [(x, y)]


    def check_collision(px, py):
        """A helper function to check if (px, py) collides with existing points
        """
        grid_px = int(px / cell_size)
        grid_py = int(py / cell_size)

        ## only search the grids around (px, py)
        for i in range(max(grid_px-2, 0), min(grid_px+3, GX)):
            for j in range(max(grid_py-2, 0), min(grid_py+3, GY)):
                if grid[i, j] > -1:
                    sx, sy = sample_list[grid[i, j]]
                    if norm2(sx-px, sy-py) < radius*radius:
                        return True
        
        return False


    ## start sampling
    while len(active_list) > 0:
        ## randomly select a point from active_list
        idx = np.random.randint(0, len(active_list))
        x, y = active_list[idx]

        success = False

        for _ in range(max_iter):
            ## sample a candidate point around (x, y)
            theta = np.random.rand() * 2 * np.pi
            dist  = (np.random.rand() + 1) * radius

            candidate_x = int(x + dist * np.cos(theta))
            candidate_y = int(y + dist * np.sin(theta))

            ## validation of candidate
            if 0 <= candidate_x < H and 0 <= candidate_y < W and grid[int(candidate_x / cell_size)][int(candidate_y / cell_size)] == -1:
                ## check collision
                if not check_collision(candidate_x, candidate_y):
                    success = True
                    break
            
        if success:
            count += 1
            sample_list.append((candidate_x, candidate_y))
            active_list.append((candidate_x, candidate_y))

            grid[int(candidate_x / cell_size), int(candidate_y / cell_size)] = count
        else:
            active_list.pop(idx)


    ## add sampled points to image
    img = np.zeros(size, np.uint8)
    for (x, y) in sample_list:
        img[x, y] = 255

    return img


def SampleDiscrepancy(img):
    """Calculate the discrepancy of a 2D image.
    """

    from scipy.stats import qmc

    loc_x, loc_y = np.where(img > 0)
    sample = np.stack([loc_x, loc_y]).T

    l_bounds = (0, 0)
    u_bounds = img.shape
    space = qmc.scale(sample, l_bounds, u_bounds, reverse=True)

    return qmc.discrepancy(space)