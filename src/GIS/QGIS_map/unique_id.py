from PyQt5.QtCore import QVariant
import pandas as pd
p=os.path.dirname(QgsProject.instance().fileName())
print(os.path.join(p,'../../../data/interm_data/mzedi_entry_skip_data_stats.csv'))

#df_mzedi = pd.read_csv('../../../data/interm_data/mzedi_entry_skip_data_stats.csv')
