#%% 
from pathlib import Path
import pandas as pd
import nibabel as nib
from nibabel.freesurfer.io import read_geometry
import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt

#%% 
project_path = Path(__file__).parent.parent

pre_implant_mri = project_path / 'data' / 'sub-RID0031' / 'derivatives' / 'freesurfer' / 'mri' / 'T1.nii.gz'
post_implant_ct = project_path / 'data' / 'sub-RID0031' / 'ses-clinical01' / 'ct' / 'sub-RID0031_ses-clinical01_acq-3D_space-T01ct_ct.nii.gz'
ct_electrodes = project_path / 'data' / 'sub-RID0031' / 'ses-clinical01' / 'ieeg' / 'sub-RID0031_ses-clinical01_space-T01ct_desc-vox_electrodes.txt'
output_dir = project_path / 'data' / 'output' / 'sub-RID0031'
freesurfer_dir = project_path / 'data' / 'sub-RID0031' / 'derivatives' / 'freesurfer'
   

file_locations_module3 = output_dir / 'ieeg_recon' / 'module3' / 'electrodes2ROI.csv'
mni305  = Path('/Users/nishant/Dropbox/Sinha/Lab/Research/t3_freesurfer/fsaverage')
mni152 = Path('/Users/nishant/Dropbox/Sinha/Lab/Research/t3_freesurfer/cvs_avg35_inMNI152')

electrodes2ROI = pd.read_csv(file_locations_module3)

talXFM_path = freesurfer_dir /  'mri' / 'transforms' / 'talairach.xfm'
t1mgz_path = freesurfer_dir / 'mri' / 'T1.mgz'

def read_talairach_xfm(fname):
    """Read the transformation matrix from a FreeSurfer .xfm file."""
    # Skip header lines until we find 'Linear_Transform'
    xfm = []
    with open(fname) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'Linear_Transform' in line:
                # Read the next 3 lines as the transformation matrix
                for j in range(3):
                    numbers = [float(x) for x in lines[i + 1 + j].strip('\n;').split()]
                    xfm.append(numbers)
    
    # Add the last row [0, 0, 0, 1] to make it a 4x4 matrix
    xfm.append([0.0, 0.0, 0.0, 1.0])
    return np.array(xfm)

# From: https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
talXFM = read_talairach_xfm(talXFM_path)
t1mgz = nib.load(t1mgz_path)
Norig = t1mgz.header.get_vox2ras()
Torig = t1mgz.header.get_vox2ras_tkr()
xform_mni305 = np.dot(talXFM, np.dot(Norig, np.linalg.inv(Torig)))
mni305_t1mgz_path = mni305 / 'mri' / 'T1.mgz'
mni305_t1mgz = nib.load(mni305_t1mgz_path)
xform_mni305_tk_ras = mni305_t1mgz.header.get_vox2ras_tkr()
xform_mni152 = np.array([[0.9975, -0.0073, 0.0176, -0.0429],
                         [0.0146, 1.0009, -0.0024, 1.5496],
                         [-0.0130, -0.0093, 0.9971, 1.1840],
                         [0.0, 0.0, 0.0, 1.0]])
mni152_t1mgz_path = mni152 / 'mri' / 'T1.mgz'
mni152_t1mgz = nib.load(mni152_t1mgz_path)
xform_mni152_tk_ras = mni152_t1mgz.header.get_vox2ras_tkr()


surfmm = electrodes2ROI.filter(['surfmm_x', 'surfmm_y', 'surfmm_z']).to_numpy()

surfmm_homog = np.hstack((surfmm, np.ones((surfmm.shape[0], 1))))
mni305_surfmm = np.round(np.dot(xform_mni305, surfmm_homog.T).T[:, :3], decimals=4)
mni305_vox = np.dot(np.linalg.inv(xform_mni305_tk_ras), 
                    np.hstack((mni305_surfmm, np.ones((mni305_surfmm.shape[0], 1)))).T).T[:, :3].astype(int)
mni305_mm = nib.affines.apply_affine(mni305_t1mgz.affine, mni305_vox)


mni305_surfmm_homog = np.hstack((mni305_surfmm, np.ones((mni305_surfmm.shape[0], 1))))

mni152_surfmm = np.round(np.dot(xform_mni152, mni305_surfmm_homog.T).T[:, :3], decimals=4)
mni152_vox = np.dot(np.linalg.inv(xform_mni152_tk_ras), 
                    np.hstack((mni152_surfmm, np.ones((mni152_surfmm.shape[0], 1)))).T).T[:, :3].astype(int)
mni152_mm = nib.affines.apply_affine(mni152_t1mgz.affine, mni152_vox)

# Create DataFrames with renamed columns for both MNI305 and MNI152
electrodes2ROI_mni305 = pd.DataFrame({
    'labels': electrodes2ROI['labels'],
    'mni305_mm_x': mni305_mm[:, 0],
    'mni305_mm_y': mni305_mm[:, 1],
    'mni305_mm_z': mni305_mm[:, 2],
    'mni305_surfmm_x': mni305_surfmm[:, 0],
    'mni305_surfmm_y': mni305_surfmm[:, 1],
    'mni305_surfmm_z': mni305_surfmm[:, 2],
    'mni305_vox_x': mni305_vox[:, 0],
    'mni305_vox_y': mni305_vox[:, 1],
    'mni305_vox_z': mni305_vox[:, 2],
    'roi': electrodes2ROI['roi'],
    'roiNum': electrodes2ROI['roiNum']
})

electrodes2ROI_mni152 = pd.DataFrame({
    'labels': electrodes2ROI['labels'],
    'mni152_mm_x': mni152_mm[:, 0],
    'mni152_mm_y': mni152_mm[:, 1],
    'mni152_mm_z': mni152_mm[:, 2],
    'mni152_surfmm_x': mni152_surfmm[:, 0],
    'mni152_surfmm_y': mni152_surfmm[:, 1],
    'mni152_surfmm_z': mni152_surfmm[:, 2],
    'mni152_vox_x': mni152_vox[:, 0],
    'mni152_vox_y': mni152_vox[:, 1],
    'mni152_vox_z': mni152_vox[:, 2],
    'roi': electrodes2ROI['roi'],
    'roiNum': electrodes2ROI['roiNum']
})

#%% 

def plot_electrodes_on_brain(coords, freesurfer_dir):
    """
    Create an interactive 3D visualization of brain surface with electrode positions
    
    Parameters:
    -----------
    electrodes2ROI : pandas.DataFrame
        DataFrame containing electrode information with columns:
        surfmm_x, surfmm_y, surfmm_z, labels, roi
    freesurfer_dir : str or Path
        Path to FreeSurfer directory containing surface files
    
    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive 3D visualization
    """
    # Convert freesurfer_dir to Path object if it's not already
    freesurfer_dir = Path(freesurfer_dir)
    
    # Load FreeSurfer surfaces
    lh_pial_verts, lh_pial_faces = read_geometry(freesurfer_dir / 'surf/lh.pial')
    rh_pial_verts, rh_pial_faces = read_geometry(freesurfer_dir / 'surf/rh.pial')

    # Create interactive 3D visualization
    fig = go.Figure()

    # Add left hemisphere mesh
    fig.add_trace(go.Mesh3d(
        x=lh_pial_verts[:, 0], y=lh_pial_verts[:, 1], z=lh_pial_verts[:, 2],
        i=lh_pial_faces[:, 0], j=lh_pial_faces[:, 1], k=lh_pial_faces[:, 2],
        color='#808080', opacity=1.0, flatshading=False,
        lighting=dict(
            ambient=0.3, diffuse=0.8, specular=1.0,
            roughness=0.1, fresnel=0.9
        ),
        lightposition=dict(x=100, y=200, z=150),
        name='Left Hemisphere',
        visible=True,
        showscale=False
    ))

    # Add right hemisphere mesh
    fig.add_trace(go.Mesh3d(
        x=rh_pial_verts[:, 0], y=rh_pial_verts[:, 1], z=rh_pial_verts[:, 2],
        i=rh_pial_faces[:, 0], j=rh_pial_faces[:, 1], k=rh_pial_faces[:, 2],
        color='#808080', opacity=1.0, flatshading=False,
        lighting=dict(
            ambient=0.3, diffuse=0.8, specular=1.0,
            roughness=0.1, fresnel=0.9
        ),
        lightposition=dict(x=100, y=200, z=150),
        name='Right Hemisphere',
        visible=True,
        showscale=False
    ))

    # Add electrode points
    fig.add_trace(go.Scatter3d(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        mode='markers',
        marker=dict(
            size=4,
            color='red',
            opacity=1.0,
            symbol='circle'
        ),
        text=electrodes2ROI['labels'],
        hovertemplate=(
            "<b>%{text}</b><br>" +
            "ROI: %{customdata[0]}<br>" +
            "X: %{customdata[1]:.2f}<br>" +
            "Y: %{customdata[2]:.2f}<br>" +
            "Z: %{customdata[3]:.2f}<br>" +
            "<extra></extra>"
        ),
        customdata=np.column_stack((
            electrodes2ROI['roi'],
            electrodes2ROI['surfmm_x'],
            electrodes2ROI['surfmm_y'],
            electrodes2ROI['surfmm_z']
        )),
        name='Electrodes',
        showlegend=True
    ))

    # Add visibility toggle buttons
    updatemenus = [
        dict(
            type="buttons",
            showactive=True,
            buttons=[
                dict(label="Show All",
                     method="update",
                     args=[{"visible": [True, True, True]}]),
                dict(label="Left Only",
                     method="update",
                     args=[{"visible": [True, False, True]}]),
                dict(label="Right Only",
                     method="update",
                     args=[{"visible": [False, True, True]}]),
                dict(label="Electrodes Only",
                     method="update",
                     args=[{"visible": [False, False, True]}]),
            ],
            direction="down",
            pad={"r": 10, "t": 10},
            x=0.02,
            xanchor="left",
            y=0.9,
            yanchor="top"
        ),
    ]

    # Update layout
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False, showgrid=False, showbackground=False),
            yaxis=dict(visible=False, showgrid=False, showbackground=False),
            zaxis=dict(visible=False, showgrid=False, showbackground=False),
            aspectmode="data",
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.5, y=1.5, z=1.5)
            ),
            bgcolor='white'
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        margin=dict(r=0, l=0, b=0, t=0),
        showlegend=False,
        updatemenus=updatemenus
    )
    
    return fig

#%% 
fig_native = plot_electrodes_on_brain(surfmm, freesurfer_dir)
fig_mni305 = plot_electrodes_on_brain(mni305_surfmm, mni305)
fig_mni152 = plot_electrodes_on_brain(mni152_surfmm, mni152)

fig_native.write_html(str(output_dir / 'ieeg_recon' / 'electrode_visualization_native.html'))
fig_mni305.write_html(str(output_dir / 'ieeg_recon' / 'electrode_visualization_mni305.html'))
fig_mni152.write_html(str(output_dir / 'ieeg_recon' / 'electrode_visualization_mni152.html'))

# %%

