
import numpy as np
import matplotlib.pyplot as plt 
import sarracen
import inspect

sdf = sarracen.read_phantom('turbbox_00067')

#info = sdf.describe()

# Calculate extra quantities 
sdf.calc_density()
sdf['vmag'] = np.sqrt(sdf['vx']**2 + sdf['vy']**2 + sdf['vz']**2)
sdf['P'] = sdf['u'] * sdf['rho'] * (sdf.params['gamma'] - 1.0)

print(sdf.columns.tolist())




fig, (ax1,ax2) = plt.subplots(ncols=2, figsize=(12,5))

ax1 = sdf.render('rho', ax=ax1, xlim=(-1, 1), ylim=(-1, 1), log_scale=True, cmap='gist_heat', vmin=1e-1, vmax=1e4, cbar_kws=dict(label='log col-dens',orientation='vertical',shrink=0.9,pad=0.03))

ax1.set_xlabel('x [pc]')
ax1.set_ylabel('y [pc]')


ax2 = sdf.render('u', ax=ax2, xlim=(-1, 1), ylim=(-1, 1), log_scale=True, cmap='plasma', vmin=1e-1, vmax=1e10, cbar_kws=dict(label='log col-etherm',orientation='vertical',shrink=0.9,pad=0.03), dens_weight=True)

ax2.set_xlabel('x [pc]')
ax2.set_ylabel(' ')


# Note: cbar_kws=dict(  --- kwargs in plt.colorbar.Colorbar()
# https://matplotlib.org/stable/api/colorbar_api.html#matplotlib.colorbar.Colorbar

'''
print(inspect.signature(sdf.render))
(target: str, x: str = None, y: str = None, z: str = None, xsec: float = None, kernel: sarracen.kernels.base_kernel.BaseKernel = None, x_pixels: int = None, y_pixels: int = None, xlim: Tuple[float, float] = None, ylim: Tuple[float, float] = None, cmap: Union[str, matplotlib.colors.Colormap] = 'gist_heat', cbar: bool = True, cbar_kws: dict = {}, cbar_ax: matplotlib.axes._axes.Axes = None, ax: matplotlib.axes._axes.Axes = None, exact: bool = None, backend: str = None, integral_samples: int = 1000, rotation: Union[numpy.ndarray, list, scipy.spatial.transform._rotation.Rotation] = None, rot_origin: Union[numpy.ndarray, list, str] = None, log_scale: bool = None, dens_weight: bool = None, normalize: bool = False, hmin: bool = False, **kwargs) -> matplotlib.axes._axes.Axes
'''

plt.show()
