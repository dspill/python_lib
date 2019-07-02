from matplotlib.widgets import RectangleSelector

 # value given by Latex \the\columnwidth
columnwidths = {
    'beamer' : 307.28987
    }

def get_figsize(columnwidth, wf=0.5, hf=(5.**0.5-1.0)/2.0, ):
    """Parameters:
      - wf [float]:  width fraction in columnwidth units
      - hf [float]:  height fraction in columnwidth units.
                     Set by default to golden ratio.
      - columnwidth [float]: width of the column in latex. Get this from LaTeX
                             using \showthe\columnwidth
    Returns:  (fig_width,fig_height): that should be given to matplotlib
    """
    fig_width_pt = columnwidth*wf
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*hf      # height in inches
    return (fig_width, fig_height)

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print('x = %d, y = %d'%(ix, iy))

    global coords
    coords.append((ix, iy))

    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)

    return coords

def fit_peak(x, y, pbl):
    max_idx = np.argmax(y)
    x0 = x[max_idx]
    y0 = y[max_idx]

    a = (x[max_idx - 1] - x0) / (y[max_idx - 1] - y0)**2

    p0 = [a, x0, y0]
    popt, pcov = curve_fit(pbl, x, y, p0=p0)

    return popt, pcov

def parabola(x, a, x0, y0):
    return a * (x - x0)**2 + y0

def toggle_selector(event):
    print('Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print('RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print('RectangleSelector activated.')
        toggle_selector.RS.set_active(True)
