import numpy as np

def microscope_objective_properties(
    make, mag, na, immersion, fn_mm=20, lambda_um=0.532, verbose=True):
    make_to_tube_lens_f_mm = {
        'Leica':200, 'Mitutoyo':200, 'Nikon':200, 'Olympus':180, 'Zeiss':165}
    tube_lens_f_mm = make_to_tube_lens_f_mm[make]
    immersion_to_n = {'air':1, 'wat':1.33, 'sil':1.41, 'gly':1.47, 'oil':1.51}
    n = immersion_to_n[immersion]
    theta = np.arcsin(na / n)
    f_mm = tube_lens_f_mm / mag
    bfp_mm = 2 * f_mm * na
    fov_um = 1000 * fn_mm / mag
    collection = 100 * (1 - np.cos(theta))  # % of hemisphere
    r0_nm = 1e3 * 0.61 * lambda_um / na     # Rayleigh resolution
    z0_nm = 1e3 * 2 * n * lambda_um / na**2 # Rayleigh axial
    xy_px_nm = r0_nm / 2                    # Nqyuist xy
    z_px_nm = z0_nm / 2                     # Nqyuist z
    Npx = 1000 * fov_um / xy_px_nm          # Number of pixels
    properties = (n, theta, f_mm, bfp_mm, fov_um, collection,
                  r0_nm, z0_nm, xy_px_nm, z_px_nm, Npx)
    if verbose:
        print('%s %ix%0.2f%s (n=%0.2f, lambda_um=%0.3f):'%(
            make, mag, na, immersion, n, lambda_um))
        print('- half angle:    %0.2fdeg (collection=%0.2f%%)'%(
            np.rad2deg(theta), collection))
        print('- focal length:  %0.2fmm (BFP=%0.2fmm)'%(f_mm, bfp_mm))
        print('- field of view: %0.2fum (FN=%0.2fmm)'%(fov_um, fn_mm))
        print('- resolution:    r0=%0.1fnm, z0=%0.1fnm (ratio: %0.1f)'%(
            r0_nm, z0_nm, z0_nm / r0_nm))
        print('- Nyquist px:    %0.1fnm(xy), %0.1fnm(z)'%(xy_px_nm, z_px_nm))
        print('- Optical px:    %i\n'%Npx)
    return properties

if __name__ == "__main__":
    properties = microscope_objective_properties('Nikon', 100, 1.35, 'sil')
