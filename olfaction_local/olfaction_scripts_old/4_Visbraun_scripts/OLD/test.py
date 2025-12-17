"""Colorer des ROI en fonction de données qui leur sont attribuées."""
import numpy as np

from visbrain.objects import RoiObj
from visbrain.utils import array2colormap

"""Tu définis un ROI object avec tes AAL."""
r_obj = RoiObj('aal')

"""Tu récupères toutes les ROI de cet atlas. C'est un pandas DataFrame.

http://visbrain.org/objects.html#visbrain.objects.RoiObj.get_labels
"""
labels = r_obj.get_labels()
# print(labels)

"""Disons que les données que tu veux associer aux ROI sont contenues dans
`data_roi`. Pour l'exmple, je ne prends que les 5 premières ROI."""
n_roi = 5  # len(labels) pour toutes les ROI (ça va être long !)
data_roi = np.arange(n_roi)  # je créé des données que l'on va associer aux ROI
labels = [38,85,1,2,34]
# print(labels)

"""Ensuite on converti ces données en couleur en utilisant une colormap. Pour
cela, il y a une fonction dans visbrain qui s'appelle `array2colormap`.

Pour plus d'info, check that (ligne 119):
https://github.com/EtienneCmb/visbrain/blob/master/visbrain/utils/color.py
"""
cmap = 'viridis'
roi_color = array2colormap(data_roi, cmap=cmap)

"""Maintenant on constuit le dictionnaire de couleurs."""
dico_color = {k: i for k, i in zip(labels, roi_color)}
# print(dico_color)

"""On sélectionne les ROI avec les couleurs qui suivent la colormap."""
r_obj.select_roi(select=labels, roi_to_color=dico_color, smooth=7)

r_obj.preview()
