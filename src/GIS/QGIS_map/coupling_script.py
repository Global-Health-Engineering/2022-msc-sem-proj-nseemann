from PyQt5.QtCore import QVariant
import pandas as pd

import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


layer=QgsVectorLayer("public-waste-skips-blantyre-malawi.shp")
field_names = layer.fields().names()
features = layer.getFeatures()
layer_provider = layer.dataProvider()
#layer_provider.addAttributes([QgsField("arrivals",QVariant.Int)])
#layer.updateFields()

df_mzedi = pd.read_csv('../../../data/interm_data/mzedi_entry_skip_data_stats.csv')
df_mzedi.set_index(df_mzedi.columns[0], inplace=True)
#print(df_mzedi.loc['Queens Guardian shelter'])
#UPDATING/ADD ATTRIBUTE VALUE
layer.startEditing()
for f in features:
    id=f.id()
    name_skip=f.attributes()[0]
    print(name_skip)
    #length=f.geometry().length()
    try:
        value_arrivals = int(df_mzedi.loc[name_skip]['#'])
    except:
        print('passed')
        value_arrivals = ''
    attr_value={6:value_arrivals}
    print(attr_value)
    layer_provider.changeAttributeValues({id:attr_value})
layer.commitChanges()