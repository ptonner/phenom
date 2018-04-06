from design import Design
import patsy

class Formula(Design):

    def __init__(self, form, name='formula', *args, **kwargs):
        super(Formula, self).__init__(name, *args, **kwargs)
        self.formula = form
        self.d = None #patsy.dmatrix(self.form, self.meta)
        self.modelDescription = patsy.ModelDesc.from_formula(self.formula)
        
    def _k(self):
        return len(self.modelDescription.rhs_termlist)

    def _matrix(self, meta):
        return np.array(patsy.dmatrix(self.modelDescription.describe(), meta))

    def _names(self):
        # return self.d.design_info.column_names
        cols = copy(self.d.design_info.column_names)

        # convert all categorical factors to something pretty
        pat = 'C\((?P<factor>[a-zA-Z0-9]+)(?P<ignore>, Treatment\([a-zA-Z0-9.]+\))?\)\[T?\.?(?P<level>[a-zA-Z0-9. ]+)\]'
        comp = re.compile(pat)
        found = map(comp.findall, cols)

        for i in range(len(cols)):
            if len(found[i]) > 0:
                cols[i] = ', '.join(['%s=%s' % (f, l) for f, _, l in found[i]])

        return cols

    def _priors(self):
        priors = -1 * np.ones(self.k)

        for t in self.d.design_info.term_names:
            priors[self.d.design_info.term_name_slices[t]] = max(priors) + 1

        return priors

    @property
    def initialized(self):
        return not self.d is None
