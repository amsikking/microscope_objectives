import numpy as np
from rayleigh_range import rayleigh_range

def intensity(r_um, z_m, zr_m, w0_um):
    r_m = 1e-6 * r_um
    w0_m = 1e-6 * w0_um
    wz_m = w0_m * (1 + (z_m / zr_m)**2)**0.5
    Irz = (w0_m / wz_m)**2 * np.exp(-2 * (r_m / wz_m)**2)
    if isinstance(r_m, int) or isinstance(r_m, float):
        print('Intensity:')
        print('(w0_um=%0.2f, zr_m=%0.2f, z_m=%0.2f, r_m=%0.2f)'%(
            w0_um, zr_m, z_m, r_m))
        print('Irz=%0.2f'%Irz)
    return Irz

if __name__ == "__main__":
    # Input:
    w0_um, n, lambda_um = 500, 1.0, 0.532
    # Get value:
    zr_m = rayleigh_range(w0_um, n, lambda_um)
    Ir0 = intensity(w0_um, 0, zr_m, w0_um)

    # Plot:
    r_m = np.linspace(-2 * w0_um, 2 * w0_um, 1000)
    Ir0 = intensity(r_m, 0, zr_m, w0_um)
    Irz = intensity(r_m, zr_m, zr_m, w0_um)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_title('Gaussian beam\n (w0_um=%0.2f, n=%0.2f, lambda_um=%0.3f)'%(
        (w0_um, n, lambda_um)))
    ax.set_ylabel('Intensity (normalised)')
    ax.set_xlabel('r (um)')
    ax.plot(r_m, Ir0, label='z=0', linestyle='-', color='b')
    ax.plot(r_m, Irz, label='z=ZR', linestyle='--', color='b')
    ax.axvline(x=w0_um, label='w0_um', linestyle='-.', color='r')
    ax.axvline(
        x=w0_um * 2**0.5, label='sqrt(2) w0_um', linestyle=':', color='r')
    plt.legend(loc="upper left")
    fig.savefig('radial_intensity.png', dpi=150)
    plt.show()
