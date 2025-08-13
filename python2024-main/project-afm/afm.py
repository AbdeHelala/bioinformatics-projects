from argparse import ArgumentParser
import pickle
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

@st.cache_data
def load_data(prefix):
    print(f"- Loading {prefix}{{.data.pickled,.heights.npy}}")
    hname = f"{prefix}.heights.npy"
    H = np.load(hname)
    m, n = H.shape
    fname = f"{prefix}.data.pickled"
    with open(fname, 'rb') as f:
        data = pickle.load(f)

    # Estimate all slopes
    slope_est = dict()  # maps tuple (s, i, j) to tuple (slope, anchors, info)
    for point, curve in data.items():
        s, i, j = point
        slope_est[point] = estimate_slope(curve, s)
    nseries = 1 + max(s for s, i, j in data.keys())
    slope_heatmaps = []
    for s in range(nseries):
        M = np.array([[(slope_est[s, i, j][0]) for j in range(n)] for i in range(m)], dtype=np.float64)
        slope_heatmaps.append(M)
    return H, data, slope_est, slope_heatmaps

def do_plot(point, curve, slope=None, anchors=None):
    """plot one distance-force curve with estimated slope"""
    s, i, j = point
    d, f = curve
    fig = plt.figure(figsize=[10, 6])
    plt.xlabel("distance (m)")
    plt.ylabel("force (N)")
    mode = 'push' if s == 0 else 'retract'
    plt.title(f"{mode} at ({i}, {j});  number of records: {len(d)}")
    label = f'data: {mode} at {(i, j)}'
    plt.scatter(d, f, s=1, label=label)
    plt.grid()
    if slope is not None and anchors is not None:
        anchor0, anchor1 = anchors[0], anchors[1]
        plt.axline(anchor0, slope=slope, color='red', linestyle='--', label=f'{slope:.4g} N/m')
        plt.plot([anchor0[0]], [anchor0[1]], 'rx')
        plt.plot([anchor1[0]], [anchor1[1]], 'rx')
    plt.legend()
    return fig

def estimate_slope(curve, s, nan=float("nan")):
    d, f = curve
    if s == 0:
        d, f = d[::-1], f[::-1]  # reverse d and f for series-0 spectra !
    # d is now increasing, f (on the left side) is decreasing.
    

    
    d_smooth = gaussian_filter1d(d, sigma=2)
    f_smooth = gaussian_filter1d(f, sigma=2)

    
    gradient = np.gradient(f_smooth, d_smooth)

    
    declining_points = gradient < 0

    
    segments = []
    segment_start = None
    for i, is_declining in enumerate(declining_points):
        if is_declining and segment_start is None:
            segment_start = i
        elif not is_declining and segment_start is not None:
            segments.append((segment_start, i - 1))
            segment_start = None
    if segment_start is not None:
        segments.append((segment_start, len(declining_points) - 1))

    
    valid_segment = max(segments, key=lambda seg: seg[1] - seg[0])

    d_segment = d_smooth[valid_segment[0]:valid_segment[1] + 1]
    f_segment = f_smooth[valid_segment[0]:valid_segment[1] + 1]

    
    if len(d_segment) < 5:
        return nan, None, None

    A = np.vstack([d_segment, np.ones(len(d_segment))]).T
    slope, intercept = np.linalg.lstsq(A, f_segment, rcond=None)[0]

    
    midpoint = len(d_segment) // 2
    anchor_distance = max(1, len(d_segment) // 6)  
    f_anchor1 = slope * d_segment[midpoint - anchor_distance] + intercept
    f_anchor2 = slope * d_segment[midpoint + anchor_distance] + intercept

    anchor1 = (d_segment[midpoint - anchor_distance], f_anchor1)
    anchor2 = (d_segment[midpoint + anchor_distance], f_anchor2)
    anchors = (anchor1, anchor2)
    info = None  # can by anything you want to return in addition
    return slope, anchors, info


# MAIN script
p = ArgumentParser()
p.add_argument("prefix",
               help="common path prefix for spectra (.data.pickled) and heights (.heights.npy)")
args = p.parse_args()
prefix = args.prefix

st.sidebar.title("AFM Data Explorer")
st.sidebar.write(f"Path prefix:\n'{prefix}'")
H, S, slope_est, slope_heatmaps = load_data(prefix)  # cached
m, n = H.shape
nseries = len(slope_heatmaps)

# Streamlit interface
series = st.sidebar.selectbox("Select series (0 for push, 1 for retract):", list(range(nseries)))
i = st.sidebar.slider("Select vertical coordinate i:", 0, m-1, 0)
j = st.sidebar.slider("Select horizontal coordinate j:", 0, n-1, 0)

st.header("Height Information")
st.write("Heatmap of height measurements.")
plt.figure(figsize=(10, 6))
plt.imshow(H, cmap='turbo', aspect='auto')
plt.colorbar(label='Height')
st.pyplot(plt)

st.header("Slope Information")
st.write(f"Heatmap of slope measurements for series {series}.")
plt.figure(figsize=(10, 6))
plt.imshow(slope_heatmaps[series], cmap='turbo', aspect='auto')
plt.colorbar(label='Slope (N/m)')
st.pyplot(plt)

point = (series, i, j)
curve = S[point]
slope, anchors, info = slope_est[point]

st.header("Measurement Series")
st.write(f"Distance vs Force curve for series {series} at point ({i}, {j}).")
fig = do_plot(point, curve, slope, anchors)
st.pyplot(fig)
