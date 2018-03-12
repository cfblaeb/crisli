import pandas as pd
from matplotlib import pyplot as plt
from sqlite3 import connect

dbc = connect("results.db")
df = pd.read_sql("select r.id, count(g.id) from region r left join grna g on r.id = g.region group by r.id", dbc)
print(df.iloc[:,1].value_counts().sort_index())
ax = df.iloc[:,1].value_counts().sort_index().plot(kind='bar', logy=True, title="gRNA per region")
ax.set_xlabel("#gRNA")
ax.set_ylabel("#regions")
plt.show()

df2 = pd.read_sql("select o.offscore from offscores o inner join grna g on o.grna = g.id", dbc)
ax = df2.plot(kind='hist', bins=100, title="offtarget score distribution", logy=True)
ax.set_xlabel("crispy++ offtarget score")
plt.show()
