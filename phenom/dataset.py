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

    def __init__(self,data, meta=None, timeColumn = 0):

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

        self.meta.index = range(self.meta.shape[0])
        self.data.columns = self.meta.index
        self.data.index.name='time'

    def concat(self, other):
        data = pd.merge(
                    self.data, other.data, left_index=True, right_index=True, how='outer')
        meta = pd.concat((self.meta, other.meta), 0)

        return DataSet(data, meta)

    def floor(self, tol=1e-9):
        self.data.values[self.data<=tol] = tol

    def copy(self):
        return DataSet(self.data, self.meta)

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

    def melt(self, norm=False, convertNames=False):

        pivot = pd.concat((self.meta, self.data.T),1,ignore_index=False)

        idvars = pivot.columns[:self.meta.shape[1]].tolist()
        valuevars = pivot.columns[self.meta.shape[1]:].tolist()
        melt = pd.melt(pivot, idvars, valuevars, var_name='time', value_name='od')
        melt.time = melt.time.astype(float)

        if norm:
            melt.od = (melt.od-melt.od.mean())/melt.od.std()

        if convertNames:
            terms = '-+'
            for c in melt.columns:
                if any([t in c for t in terms]):
                    temp = c
                    for t in terms:
                        temp = temp.replace(t, '_')
                    melt[temp] = melt[c]
                    del melt[c]

        return melt

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

    def poly_scale(self,p,ind=5,groupby=None):
        """Scale growth data by a polynomial of degree p, using the first ind datapoints, grouping by groupby."""

        time = self.data.index.values[:ind]

        if groupby is None:
            #group = {(None,self.meta.index)}
            for i in range(self.data.shape[1]):
                temp = self.data.iloc[:,i]
                od = temp.values[:ind]

                coeff = np.polyfit(time,od,p)

                temp = temp - np.polyval(coeff,self.data.index.values[0])
                self.data.iloc[:,i] = temp
        else:
            group = self.meta.groupby(groupby)
            for k,index in group.groups.items():
                temp = self.data.loc[:,index]
                od = temp.values[:ind,:].ravel()

                coeff = np.polyfit(time.tolist()*temp.shape[1],od,p)

                temp = temp - np.polyval(coeff,self.data.index.values[0])
                self.data.loc[:,index] = temp

    def filter(self):
        """Remove data rows where observations are missing (in any column!)"""
        self.data = self.data.loc[~self.data.isnull().any(1),:]

    def plot(self, columns=[], rows=[], colorby=[], buff=.1, colorLabels=True):

        import matplotlib.pyplot as plt

        singleTerm = False
        if len(columns) + len(rows) + len(colorby) == 1:
            singleTerm = True

        cmap = plt.get_cmap()
        labels = set()
        g = self.meta.groupby(columns + rows + colorby)

        # legend variables
        handles = []
        lines = []

        axes = {}

        columnKeys = rowKeys = colorKeys = []
        if len(columns) > 0:
            columnKeys = list(self.meta.groupby(columns).groups.keys())
            ncol = len(columnKeys)
        else:
            ncol = 1

        if len(rows) > 0:
            rowKeys = list(self.meta.groupby(rows).groups.keys())
            nrow = len(rowKeys)
        else:
            nrow = 1

        if len(colorby) > 0:
            colorKeys = list(self.meta.groupby(colorby).groups.keys())


        keys = g.groups.keys()
        keys = sorted(keys)
        for k in keys:
        # for k, ind in g.groups.iteritems():

            ind = g.get_group(k).index

            cval = k[:len(columns)]
            rval = k[len(columns):len(columns)+len(rows)]
            colorval = k[len(columns)+len(rows):]

            # this is only okay if there is more than one column contributing to keys
            if not singleTerm:
                if len(columns) == 1:
                    cval = cval[0]
                if len(rows) == 1:
                    rval = rval[0]
                if len(colorby) == 1:
                    colorval = colorval[0]

            i = j = 0
            if len(columnKeys):
                i = columnKeys.index(cval)
            if len(rowKeys):
                j = rowKeys.index(rval)

            # ax = plt.subplot(nrow, ncol, j*ncol + i + 1)

            if (j, i) in axes:
                #plt.axis(axes[(i,j)])
                ax = axes[(j,i)]
            else:
                ax = plt.subplot2grid((nrow+1, ncol + 1), (j+1, i+1))
                axes[(j,i)] = ax

            title = ""
            if i == 0:

                if len(rows) > 1:
                    s = ['%s = %s'%(r, str(v)) for r,v in zip(rows, rval)]
                    title +=  ', '.join(s)
                    title += ' | '
                elif len(rows) == 1:
                    title +=  '%s' % (str(rval))
                    title += ' | '

            if len(columns) > 1:
                s = ['%s = %s'%(r, str(v)) for r,v in zip(columns, cval)]
                title +=  ', '.join(s)
            elif len(columns) == 1:
                title +=  '%s' % (str(cval))

            #plt.title(title)
            # ax.set_title(title)

            if len(colorby):
                f = 1.*(colorKeys.index(colorval) + buff) / (len(colorKeys) + 2*buff)
                c = cmap(f)

                if colorval in labels:
                    # plt.plot(self.data.index, self.data.iloc[:,ind],color=c)
                    ax.plot(self.data.index, self.data.iloc[:,ind],color=c)
                else:
                    labels.add(colorval)

                    l = ""
                    if len(colorby) == 1:
                        l += colorby[0]
                    else:
                        l += ', '.join(colorby)
                    l+= ' = '
                    if len(colorby) == 1:
                        l += str(colorval)
                    else:
                        l += ', '.join(map(str, colorval))

                    # '%s=%s'%(', '.join(colorby), ', '.join(col))

                    # plt.plot(self.data.index, self.data.iloc[:,ind[0]], color=c, label=l)
                    # plt.plot(self.data.index, self.data.iloc[:,ind[1:]],color=c)
                    ax.plot(self.data.index, self.data.iloc[:,ind[0]], color=c, label=l)
                    ax.plot(self.data.index, self.data.iloc[:,ind[1:]],color=c)
            else:
                # plt.plot(self.data.index, self.data.iloc[:,ind])
                ax.plot(self.data.index, self.data.iloc[:,ind])

            mi = self.data.min().min()
            ma = self.data.max().max()
            ran = ma - mi
            # plt.ylim(mi - .02*ran, ma + .02*ran)
            ax.set_ylim(mi - .02*ran, ma + .02*ran)

            for h, l in zip(*ax.get_legend_handles_labels()):
                if not h in handles:
                    handles += [h]
                    lines += [l]

            #if colorLabels:
            #    plt.legend()

        for i, c in enumerate(columnKeys):
            ax = plt.subplot2grid((nrow+1, ncol + 1), (0, i+1))
            ax.text(0.5,0.5,c,horizontalalignment='center')
            ax.set_axis_off()


        for i, r in enumerate(rowKeys):
            ax = plt.subplot2grid((nrow+1, ncol + 1), (i+1,0))
            ax.text(0.5,0.5,r,horizontalalignment='center')
            ax.set_axis_off()

        if colorLabels:
            # ax = plt.subplot(nrow, ncol, 1)
            ax = plt.subplot2grid((nrow+1, ncol + 1), (0, 0))
            ax.legend(handles, lines)
            ax.set_axis_off()
