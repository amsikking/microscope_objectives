import numpy as np
import scipy.special as sp

def psf_properties(na, immersion, r_um, z_um, lambda_um, verbose=True):
    immersion_to_n = {'air':1, 'wat':1.33, 'sil':1.41, 'gly':1.47, 'oil':1.51}
    n = immersion_to_n[immersion]
    v = (2 * np.pi * na / lambda_um) * r_um
    u = (2 * np.pi * na**2 / (n * lambda_um)) * z_um
    r0_um = 0.61 * lambda_um / na           # Rayleigh resolution
    z0_um = 2 * n * lambda_um / na**2       # Rayleigh axial
    with np.errstate(divide='ignore', invalid='ignore'): # avoid div 0 error
        ir = (2 * sp.j1(v) / v)**2          # normalized radial intensity
    iz = np.sinc((u/4)/np.pi)**2            # normalized axial intensity
    pwr = 1 - sp.j0(v)**2 - sp.j1(v)**2     # normalized power thru aperture
    if verbose:
        if np.isnan(ir): ir = 1 # fix div 0
        print('PSF properties ', end='')
        print('(%0.2f NA, n=%0.2f, lambda_um=%0.2f):'%(na, n, lambda_um))
        print('- r0=%0.2fnm, z0=%0.2fnm (ratio=%0.2f)'%(
            1e3 * r0_um, 1e3 * z0_um, z0_um / r0_um))
        print('- Depth of field: %0.2fnm'%(1e3 * z0_um / 2))
        print('- Radial intensity (r = %0.2fum): %0.2f%% of peak'%(
            r_um, 100 * ir))
        print('- Axial  intensity (z = %0.2fum): %0.2f%% of peak \n'%(
            z_um, 100 * iz, ))
    return (n, r0_um, z0_um, ir, iz, pwr)

if __name__ == "__main__":
    # Input:
    na, imm, r_um, z_um, lambda_um = 1.35, 'sil', 0, 0, 0.532
    # Get value:
    n, r0_um, z0_um, ir, iz, pwr = psf_properties(
        na, imm, r_um, z_um, lambda_um)

    # Plot:
    d_um = np.linspace(0, 2 * z0_um, 1000)
    n, r0_um, z0_um, ir, iz, pwr = psf_properties(
        na, imm, d_um, d_um, lambda_um, verbose=False)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_title('Point spread function\n' +
                 '(NA=%0.2f, n=%0.2f, lambda_um=%0.3f, z0/r0=%0.2f)'%(
                     (na, n, lambda_um, z0_um/r0_um)))
    ax.set_xlabel('Displacement (um)')
    ax.set_ylabel('Intensity/power (normalised)')
    ax.plot(d_um, ir, 'r', label='Radial')
    ax.axvline(x=r0_um, label='r0 (%0.3fum)'%(r0_um), color='r')    
    ax.plot(d_um, iz, 'b', label='Axial', linestyle='--')
    ax.axvline(x=z0_um, label='z0 (%0.3fum)'%(z0_um), linestyle='--', color='b')
    ax.plot(d_um, pwr, 'k', label='Power thru aperture', linestyle=':')
    ax.axvline(x=(z0_um/4), label='z tol. 80%% (%0.3fum)'%(z0_um/4),
                linestyle='-.', color='k')
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)
    ax.legend(loc = "upper right", framealpha=1 )
    fig.savefig('psf.png', dpi=150)
    plt.show()
