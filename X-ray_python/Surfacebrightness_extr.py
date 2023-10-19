import pyproffit

dat=pyproffit.Data(imglink='/Users/thomas_1/Desktop/A1446/4975/repro/flux/broad_thresh.img',
                   explink='/Users/thomas_1/Desktop/A1446/4975/repro/flux/broad_thresh.expmap')

prof=pyproffit.Profile(dat,center_choice='peak',maxrad=45.,binsize=20.,centroid_region=30.)
prof.SBprofile(ellipse_ratio=1)

mod=pyproffit.Model(pyproffit.BetaModel)
fitobj=pyproffit.Fitter(model=mod, profile=prof, beta=0.7, rc=2.,norm=-2,bkg=-4)
fitobj.Migrad()
prof.Backsub(fitobj)
prof.Plot(model=mod,outfile='surbri.pdf')
