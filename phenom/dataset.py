import pandas as pd
import numpy as np

nan_or_zero = lambda x: np.isnan(x) or abs(x) < 1e-9

class DataSet(object):

    @classmethod
    def fromDirectory(cls, dir, datafile = 'data.csv', metafile = 'meta.csv', *args, **kwargs):
        import os

        assert datafile in os.listdir(dir)
        assert metafile in os.listdir(dir)

        data = pd.read_csv(os.path.join(dir, datafile))
        meta = pd.read_csv(os.path.join(dir, metafile))

        return cls(data, meta, *args, **kwargs)

    def __init__(self,data,meta=None, timeColumn = 0):

        if not type(data) == pd.DataFrame:
            data = pd.DataFrame(data)

        if meta is None:
            meta = pd.DataFrame(index=range(data.shape[1]))

        if not type(meta) == pd.DataFrame:
            meta = pd.DataFrame(meta)

        self.meta = meta
        self.data = data

        if self.data.shape[1] == self.meta.shape[0] + 1:
            self.data.index = self.data.iloc[:,timeColumn]
            self.data = self.data.drop(self.data.columns[timeColumn],1)

        assert self.data.shape[1] == self.meta.shape[0], 'frames do no match, %d x %d' % (self.data.shape[1], self.meta.shape[0])

        self.data.columns = self.meta.index
        self.data.index.name='time'

    def copy(self):
        return DataSet(self.data.copy(), self.meta.copy())

    def __eq__(self,other):
        if not type(other) == DataSet:
            return False

        if any([x not in self.meta.columns for x in other.meta.columns]):
            return False

        if not other.data.shape == self.data.shape or not other.meta.shape == self.meta.shape:
            return False

        other = other.copy()
        other.meta = other.meta[self.meta.columns]

        # return ((self.data-other.data).applymap(nan_or_zero)).all().all() and (self.meta == other.meta).all().all()
        return self.data.equals(other.data) and self.meta.equals(other.meta)#(self.meta == other.meta).all().all()

    def __repr__(self,):
        return 'Dataset: %d x %d, %d x %d\n%s\n%s' % (self.data.shape[0], self.data.shape[1], self.meta.shape[0], self.meta.shape[1], str(self.data.head()), str(self.meta.head()))

    def trim(self, start, stop=None):
        """remove samples before start and after stop."""
        if stop is None:
            stop = self.data.shape[0]

        start = max(start, 0)
        stop = min(stop, self.data.shape[0])
        self.data = self.data.iloc[start:stop,:]

    def log(self):

        self.data[self.data<=0] = np.nan

        self.data = np.log2(self.data)

    def melt(self):
	"""melt this dataset into 'tidy' format, where each time-point of an individual replicate is a seperate row."""
        pivot = pd.concat((self.meta, self.data.T),1,ignore_index=False)

        idvars = pivot.columns[:self.meta.shape[1]].tolist()
        valuevars = pivot.columns[self.meta.shape[1]:].tolist()
        return pd.melt(pivot, idvars, valuevars, var_name='time', value_name='od')

    def build(self,effects=[],covariates=[],scale=None,**kwargs):

        if 'x' in covariates:
            covariates.remove('x')

        if len(covariates)>0:
            temp = pd.concat((self.meta,self.data.T),1)
            temp['rep'] = temp.index

            tidy = pd.melt(temp,self.meta.columns.tolist()+['rep'],
                            self.data.index.tolist(),
                            var_name='x',value_name='y')

            pivot = pd.pivot_table(tidy,values='y',
                            index=['x']+covariates,
                            columns=['rep'])

            x = pd.DataFrame(np.array(pivot.index.tolist()),columns=['x']+covariates).values
            y = pivot.values

        else:
            x = pd.DataFrame(self.data.index.values,columns=['x']).values
            y = self.data.values

        effect = self.meta[effects]
        labels = []

        select = [True]*self.meta.shape[0]
        for k in kwargs.keys():
            if k in self.meta:
                if type(kwargs[k]) == list:
                    select = (select) & (self.meta[k].apply(lambda x: x in kwargs[k]))
                else:
                    select = (select) & (self.meta[k] == kwargs[k])
        y = y[:,np.where(select)[0]]
        effect = effect.loc[select,:]

        for e in effect.columns:
            temp,l = pd.factorize(effect[e])
            effect[e] = temp
            labels.append(l.tolist())

        if scale=='range':
            x = (x-x.min())/(x.max()-x.min())
        elif scale=='norm':
            x = (x-x.mean())/x.std()

        return x,y,effect,labels
