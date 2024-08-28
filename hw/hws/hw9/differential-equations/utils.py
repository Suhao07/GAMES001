import numpy as np

def compute_laplace_residual(u):
    residual = 0
    grid_size = u.shape
    for i in range(1, grid_size[0]-1):
        for j in range(1, grid_size[1]-1):
            residual += (u[i-1, j] + u[i+1, j] + u[i, j-1] + u[i, j+1] - 4 * u[i, j]) ** 2
    return np.sqrt(residual)

def iter_jacobi(u):
    # TODO: your code here
    # raise NotImplementedError
    H, W = u.shape
    u_prev = u.copy()
    
    for i in range(1, H-1):
        for j in range(1, W-1):
            u[i, j] = 1/4 * (u_prev[i-1, j] + u_prev[i+1, j] + u_prev[i, j-1] + u_prev[i, j+1])

def explicit_euler(position, velocity, acceleration, dt):
    # TODO: your code here
    # raise NotImplementedError

    ## update position
    _position = position + velocity * dt

    ## update velocity
    _velocity = velocity + acceleration(position) * dt

    return _position, _velocity

def symplectic_euler(position, velocity, acceleration, dt):
    # TODO: your code here
    # raise NotImplementedError

    ## update velocity
    velocity = velocity + acceleration(position) * dt

    ## update position
    position = position + velocity * dt

    return position, velocity

def symplectic_euler_2(position, velocity, acceleration, dt):
    # TODO: your code here
    # raise NotImplementedError

    ## update mid-point velocity
    velocity_mid = velocity + 0.5 * acceleration(position) * dt
    
    ## update position
    position = position + velocity_mid * dt

    ## update velocity
    velocity = velocity_mid + 0.5 * acceleration(position) * dt

    return position, velocity

def symplectic_euler_3(position, velocity, acceleration, dt):
    # TODO: your code here
    # raise NotImplementedError

    def expDV(d):
        return velocity + d * acceleration(position) * dt
    
    def expDT(c):
        return position + c * velocity * dt
    
    ci = [1, -2/3, 2/3]
    di = [-1/24, 3/4, 7/24]

    ci.reverse()
    di.reverse()

    ## update velocity and position iteratively
    for c, d in zip(ci, di):
        velocity = expDV(d)
        position = expDT(c)

    return position, velocity
