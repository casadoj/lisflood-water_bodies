import numpy as np
import pandas as pd
from xml.dom import minidom
from typing import Union, List, Tuple, Dict
from pathlib import Path



def try_float(string):
    try:
        return float(string)
    except ValueError:
        return string
        
        

def reservoir_attributes(file: str, name: Union[str, int] = None) -> pd.DataFrame:
    """Given an XML file with the attributes of a reservoir in the SNCZI database (Sistema Nacional de Cartografía de Zonas Inundables), it extracts the data and organizes it in a table
    
    Inputs:
    -------
    file: str
        File name (including extension XML) to be read
    name: Union[str, int]
        Identifier of the reservoir that will be the column name in the resulting table
        
    Output:
    -------
    attributes: pd.DataFrame
        Table with a single column with the attributes of the reservoir
    """
    
    # read file and extract headers and values of the attribute table
    headers, values = [], []
    xmldoc = minidom.parse(file)
    tables = xmldoc.getElementsByTagName('Table')
    for t, table in enumerate(tables):
        rows = table.getElementsByTagName('Row')
        for r, row in enumerate(rows):
            if r == 0:
                continue
            cells = row.getElementsByTagName('Cell')
            #print(f'row {i}\t{len(cells)} cells')
            for c, cell in enumerate(cells):
                for child in cell.childNodes:
                    if r == 1:
                        headers.append(child.childNodes[0].nodeValue)
                    elif r == 2:
                        try:
                            values.append(child.childNodes[0].nodeValue)
                        except:
                            values.append('---')

    # adapt values: convert to float (if possible) and divide the coordinates in two values
    values = [v.replace('.', '').replace(',', '.') if (',' in v) or ('.' in v) else v for v in values]
    data = []
    for value in values:
        if len(value.split(' - ')) == 2:
            try:
                data += [float(x) for x in value.split(' - ')]
            except:
                data.append(value)
        else:
            data.append(try_float(value))

    # create pandas.DataFrame of the attributes
    attributes = pd.DataFrame(pd.Series(data=data, index=headers, name=name, dtype='object'))
    attributes = attributes.replace('---', np.nan)
    attributes = attributes[~attributes.index.duplicated(keep='first')]
    #idx_int = ['Código del embalse',
    #           'Coord. X UTM ETRS89 Huso 30',
    #           'Coord. Y UTM ETRS89 Huso 30',
    #           'Coord. X Manual',
    #           'Coord. Y Manual',
    #           'Id. Hoja 1:50.000',
    #           'Código de infraestructura']
    #attributes.loc[idx_int] = attributes.loc[idx_int].astype(int)

    return attributes



def dam_attributes(file: str, name: Union[str, int] = None) -> pd.DataFrame:
    """Given an XML file with the attributes of a dam in the SNCZI database (Sistema Nacional de Cartografía de Zonas Inundables), it extracts the data and organizes it in a table
    
    Inputs:
    -------
    file: str
        File name (including extension XML) to be read
    name: Union[str, int]
        Identifier of the dam that will be the column name in the resulting table
        
    Output:
    -------
    attributes: pd.DataFrame
        Table with a single column with the attributes of the dam
    """
    
    # read file and extract headers and values of the attribute table
    xmldoc = minidom.parse(file)
    headers, values = [], []
    rows = xmldoc.getElementsByTagName('Row')
    for i, row in enumerate(rows):
        cells = row.getElementsByTagName('Cell')
        if len(cells) == 1:
            continue
        else:
            for cell in cells:
                for child in cell.childNodes:
                    if i == 1:
                        headers.append(child.childNodes[0].nodeValue)
                    elif i == 2:
                        values.append(child.childNodes[0].nodeValue)

    # adapt headers so that the coordinates have two separate fields
    headers = [x.strip()[:-1] for x in headers]
    index = []
    for x in headers:
        if 'Coordenadas' not in x:
            index.append(x)
        else:
            index += ['X-UTM30ETRS89', 'Y-UTM30ETRS89']

    # adapt values: convert to float (if possible) and divide the coordinates in two values
    values = [v.replace('.', '').replace(',', '.') if (',' in v) or ('.' in v) else v for v in values]
    data = []
    for value in values:
        if len(value.split(' - ')) == 2:
            try:
                data += [float(x) for x in value.split(' - ')]
            except:
                data.append(value)
        else:
            data.append(try_float(value))

    # create pandas.Series
    attributes = pd.DataFrame(pd.Series(data=data, index=index, name=name, dtype='object'))
    attributes = attributes.replace('---', np.nan)
    
    return attributes