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
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import AnchoredText
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


class Box:
    def __init__(self, ax, point, x, y, width: float = 1.0, height: float = 1.0, color='black'):
        self.ax = ax
        self.point = point
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.color = color
    
    def draw(self):
        box = Rectangle((self.x, self.y), self.width, self.height, ec=self.color, fc='white', zorder=5)
        self.ax.add_patch(box)
        line = Line2D([self.point[0], self.x + self.width / 2], [self.point[1], self.y], lw=1, color='gray')
        self.ax.add_line(line)


def position_cards(ax, scores: list):
    """
    A helper function for position_cards near the data point as close as possible
    """
    
    left, right = ax.get_xlim()
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
    
    def __init__(self, ax, pos, width, point, smiles, data):
        self.ax = ax
        self.x = pos
        self.width = width
        self.point = point
        self.smiles = smiles
        self.data = data
        
    def draw(self):
        y = self.ax.get_ylim()[1]
        for i, smiles in enumerate(self.smiles):
            image = smiles_image(smiles, legend=i == 3)
            box = AnnotationBbox(
                OffsetImage(image, zoom=0.30),
                xy=(0.5, 0), xybox=(self.x, y * (0.4 + i * 0.15)),
                xycoords='axes fraction', boxcoords="data",
                box_alignment=(-0.02, 0.01), pad=0, frameon=False, zorder=10
            )
            self.ax.add_artist(box)
        
        
            
        Box(self.ax, self.point, self.x, y * 0.2, width=self.width, height=y * 0.795).draw()
    
    def __str__(self):
        return f'{self.data.get("compound", "Compound")} {self.pos} {self.smiles[0]}'


def overview_plot(df: pl.DataFrame, **kwargs):
    image = kwargs['outdir'].joinpath('overview.png')
    libraries = kwargs['libraries']
    num, cols = len(libraries), 6
    rows, mod = divmod(num, cols)
    rows = rows + 1 if mod else (rows or 1)
    low, high = df[kwargs['xs']].min(), df[kwargs['ys']].max()
    low, high = low * 0.9 if low >= 0 else low * 1.1, high * 1.1
    
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(kwargs['figure_width'], kwargs['figure_height']))
    for i, ax in enumerate(axes.flatten()):
        row, col = divmod(i, cols)
        
        try:
            library = libraries[i]
        except IndexError:
            axes[row - 1][col].xaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.set_xticks([]), ax.set_xticklabels([])
            ax.set_yticks([]), ax.set_yticklabels([])
            ax.get_xaxis().set_visible(False), ax.get_yaxis().set_visible(False)
            _ = [spine.set_visible(False) for spine in ax.spines.values()]
            continue
        
        dl = df.filter(pl.col('library') == library)
        xs, ys, cs = kwargs['xs'], kwargs['ys'], kwargs['circle_sizes'][0]
        
        dm = dl.filter(pl.col('axis').is_in([0, 1, 2]))
        ax.scatter(dm[xs], dm[ys], c=dm['color'], s=cs, edgecolors='#6f6f6f', lw=0.3, zorder=2, clip_on=False)
        
        dd = dl.filter(pl.col('axis').is_in([3, 4, 5]))
        ax.scatter(dd[xs], dd[ys], c=dd['color'], s=cs, edgecolors='#6f6f6f', lw=0.3, zorder=4, clip_on=False)
        
        dt = dl.filter(pl.col('axis') == 6)
        ax.scatter(dt[xs], dt[ys], c=dt['color'], s=cs, edgecolors='#6f6f6f', lw=0.3, zorder=6, clip_on=False)
        
        color_map = {'Strong': '#E60003', 'Weak': '#FF7375'}
        color = next((color_map[k] for k in color_map if k in set(dl['enrichment'])), 'black')
        ax.set_title(library, fontweight='bold', color=color, y=0.88)
        
        ax.set_xlim(low, high), ax.set_ylim(low, high)
        ax.margins(0.05)
        
        if row == rows - 1:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        else:
            ax.set_xticks([])
        if col == 0:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        else:
            ax.set_yticks([])
    
    fig.text(0.5, 0.01, kwargs['x_label'], ha='center', va='center', fontsize=12)
    fig.text(0.02, 0.5, kwargs['y_label'], ha='center', va='center', rotation='vertical', fontsize=12)
    
    fig.subplots_adjust(left=0.05, right=0.99, top=0.97, bottom=0.05, wspace=0.06, hspace=0.08)
    fig.savefig(image, dpi=kwargs['dpi'])
    logger.debug(f'Overview plot was saved to {image}\n')