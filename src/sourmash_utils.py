"Python utilities for sourmash plugins and scripts."
import sourmash

__all__ = ['FracMinHash']

class FracMinHash(sourmash.MinHash):
    def __init__(self, *,
                 ksize=31,
                 scaled=1000,
                 moltype='DNA',
                 track_abundance=False,
                 **kwargs):

        kwargs = dict(kwargs)
        assert 'is_protein' not in kwargs
        assert 'dayhoff' not in kwargs
        assert 'hp' not in kwargs

        if moltype == 'DNA':
            kwargs['is_protein'] = False
            kwargs['dayhoff'] = False
            kwargs['hp'] = False
        elif moltype == 'protein':
            kwargs['is_protein'] = True
        elif moltype == 'dayhoff':
            kwargs['dayhoff'] = True
        elif moltype == 'hp':
            kwargs['hp'] = True
        else:
            raise ValueError(f"unknown moltype: '{moltype}'")

        super().__init__(n=0,
                         ksize=ksize,
                         scaled=scaled,
                         track_abundance=track_abundance,
                         **kwargs)
