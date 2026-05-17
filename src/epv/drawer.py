import io
import sys
from pathlib import Path
from functools import partial
from datetime import datetime
from collections.abc import Sequence

from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw, Descriptors

import polars as pl
from PIL import Image
from loguru import logger
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, DrawingArea, TextArea

OPTIONS = {
    'dpi': {'help': 'The DPI of the output plot', 'default': 300, 'type': int},
    'figure_width': {'help': 'The width of figure in inch', 'default': 18, 'type': int},
    'figure_height': {'help': 'The height of figure in inch', 'default': 9, 'type': int},
    'circle_sizes': {'help': 'The circle sizes of overview and individual library scatter plots',
                     'default': [50, 500], 'type': Sequence},
    'outdir': {'help': 'Path string to a output directory'},
    'verbose': {'help': 'Enable verbose mode to emit verbose messages', 'type': bool, 'default': False}
}


def position_cards(ax, scores: list):
    """
    A helper function for position_cards near the data point as close as possible
    """
    
    # TODO: fine tune the right padding
    left, right = ax.get_xlim()
    print(left, right)
    n, xs = len(scores), []
    pad = (right - left) * 0.01
    width = (right - left - pad * 7) / 6
    half = width / 2
    
    for i, score in enumerate(scores):
        x = score - half - pad
        if score < left:
            x = left + pad
        elif score + half > right:
            x = right - width - pad
        if i == 0:
            xs.append(x + pad)
        else:
            if xs:
                previous = xs[-1]
                if previous + pad + width > x:
                    offset = previous + pad + width - x
                    if offset + x + pad + width > right:
                        xs = [p - offset for p in xs]
                        xs.append(x)
                    else:
                        xs.append(x+offset)
                else:
                    xs.append(x)
            else:
                xs.append(x)
    return xs, width
    
    
def smiles_image(smiles: str, width: int = 500, height: int = 200,
                 bond_length=15, bond_line_width=1.0, legend=False):
    mol = Chem.MolFromSmiles(smiles.split()[0])
    if mol:
        drawer = Draw.MolDraw2DCairo(500, 200)
        options = drawer.drawOptions()
        options.bondLineWidth = bond_line_width
        options.fixedBondLength = bond_length
        options.additionalAtomLabelPadding = 0.05
        options.addAtomIndices = False
        options.includeAtomTags = False
        options.addStereoAnnotation = True
        
        Chem.rdDepictor.Compute2DCoords(mol)
        Draw.rdMolDraw2D.PrepareMolForDrawing(mol)
        
        if legend:
            mw = Descriptors.MolWt(mol)
            alogp = Descriptors.MolLogP(mol)
            legend = f"MW: {mw:.0f}, ALogP: {alogp:.2f}"
            options.useBWAtomPalette()
            options.legendFraction = 0.2
        else:
            legend = ''
        drawer.DrawMolecule(mol, legend=legend)
        drawer.FinishDrawing()
        image_data = drawer.GetDrawingText()
        image = Image.open(io.BytesIO(image_data))
        drawer.ClearDrawing()
    else:
        image = np.zeros((height, width, 4), dtype=np.uint8)
        image[:, :, 3] = 0
    return image


class Card:
    """A compound card"""
    
    def __init__(self, compound, smiles, data, point, x, y, width, ax, fig):
        self.compound = compound
        self.smiles = smiles
        self.data = data
        self.point = point
        self.x = x
        self.y = y
        self.width = width
        self.ax = ax
        self.fig = fig
        
    def render(self):
        ab = AnnotationBbox(
                DrawingArea(150, 335, 0, 0),
                self.point,
                xybox=(self.x, self.y * 0.2),
                xycoords='data',
                boxcoords="data",
                box_alignment=(0, 0)
        )
        self.ax.add_artist(ab)

        invert = self.ax.transData.inverted()

        renderer = self.fig.canvas.get_renderer()
        box = ab.get_window_extent(renderer=renderer)
        bottom_center = invert.transform((box.x0 + box.width * 0.5, box.y0))
        arrow = FancyArrowPatch(self.point, bottom_center, arrowstyle="-", color='gray', lw=0.5)
        self.ax.add_patch(arrow)
        
        y = self.ax.get_ylim()[1]
        for i, smiles in enumerate(self.smiles):
            image = smiles_image(smiles, legend=i == 3)
            box = AnnotationBbox(
                OffsetImage(image, zoom=0.3),
                xy=(0.5, 0),
                xybox=(self.x, y * (0.4 + i * 0.15)),
                xycoords='axes fraction',
                boxcoords="data",
                box_alignment=(0, 0),
                pad=0, frameon=True
            )
            self.ax.add_artist(box)
    
        
    def __str__(self):
        return f'{self.compound or "Compound"} {self.pos} {self.smiles[0]}'
        