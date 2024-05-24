from pyfv3jedilm import fv3jedi_lm_mod

import torch
 
# A (Work-In-Progress) PyTorch interface to FV3 Linear Model  

class FV3LMFnc(torch.autograd.Function):
    @staticmethod
    def forward(lm, u, v, t, delp, tracers, w=None, delz=None,
                cfcn=None, qls=None, qcn=None, phis=None, frocean=None, frland=None, varflt=None, ustar=None,
                bstar=None, zpbl=None, cm=None, ct=None, cq=None, kcbl=None, ts=None, khl=None, khu=None):
        
        # differentiable trajectory args
        lm.traj.u = u.detach().numpy()
        lm.traj.v = v.detach().numpy()
        lm.traj.t = t.detach().numpy()
        lm.traj.delp = delp.detach().numpy()
        lm.traj.tracers = tracers.detach().numpy()

        
        FV3LMFnc.set_optional_args(lm, w=w, delz=delz, cfcn=cfcn, qls=qls, qcn=qcn, phis=phis, frocean=frocean,
                          frland=frland, varflt=varflt, ustar=ustar, bstar=bstar, zpbl=zpbl, cm=cm,
                          ct=ct, cq=cq, kcbl=kcbl, ts=ts, khl=khl, khu=khu)
        
        lm.step_nl()

        result = (u.new(lm.traj.u),
                  v.new(lm.traj.v),
                  t.new(lm.traj.t),
                  delp.new(lm.traj.delp),
                  tracers.new(lm.traj.tracers))
        
        if w is not None:
            result = (*result, w.new(lm.traj.w))

        if delz is not None:
            result = (*result, delz.new(lm.traj.delz))

        if cfcn is not None:
            result = (*result, cfcn.new(lm.traj.cfcn))

        return result

    @staticmethod
    def setup_context(ctx, inputs, output):
        lm = inputs[0]
        input_args = inputs[1:]

        ctx.lm = lm
        ctx.save_for_backward(*input_args)

    @staticmethod
    def set_optional_args(obj, **kwargs):
        for key, value in kwargs.items():
            if value is not None:
                setattr(obj, key, value.detach().numpy())

    @staticmethod
    def backward(ctx, Du, Dv, Dt, Ddelp, Dtracers, Dw=None, Ddelz=None, Dcfcn=None):
        
        lm = ctx.lm
        
        (u, v, t, delp, tracers, w, delz,
         cfcn, qls, qcn, phis, frocean, frland, varflt, ustar,
         bstar, zpbl, cm, ct, cq, kcbl, ts, khl, khu) = ctx.saved_tensors

        # differentiable trajectory args
        lm.traj.u = u.detach().numpy()
        lm.traj.v = v.detach().numpy()
        lm.traj.t = t.detach().numpy()
        lm.traj.delp = delp.detach().numpy()
        lm.traj.tracers = tracers.detach().numpy()

        FV3LMFnc.set_optional_args(lm, w=w, delz=delz, cfcn=cfcn, qls=qls, qcn=qcn, phis=phis, frocean=frocean,
                                   frland=frland, varflt=varflt, ustar=ustar, bstar=bstar, zpbl=zpbl, cm=cm,
                                   ct=ct, cq=cq, kcbl=kcbl, ts=ts, khl=khl, khu=khu)

        # perturbation args
        lm.pert.u = Du.numpy()
        lm.pert.v = Dv.numpy()
        lm.pert.t = Dt.numpy()
        lm.pert.delp = Ddelp.numpy()
        lm.pert.tracers = Dtracers.numpy()

        if Dw is not None:
            lm.pert.w = Dw.numpy()

        if Ddelz is not None:
            lm.pert.delz = Ddelz.numpy()

        if Dcfcn is not None:
            lm.pert.cfcn = Dcfcn.numpy()

        lm.step_ad()

        return (None, # lm
                Du.new(lm.pert.u),
                Dv.new(lm.pert.v),
                Dt.new(lm.pert.t),
                Ddelp.new(lm.pert.delp),
                Dtracers.new(lm.pert.tracers),
                Dw.new(lm.pert.w) if Dw is not None else None,
                Ddelz.new(lm.pert.delz) if Ddelz is not None else None,
                Dcfcn.new(lm.pert.cfcn) if Dcfcn is not None else None,
                None, # qls
                None, # qcn
                None, # phis
                None, # frocean
                None, # frland
                None, # varflt
                None, # ustar
                None, # bstar
                None, # zpbl
                None, # cm
                None, # ct
                None, # cq
                None, # kcbl
                None, # ts
                None, # khl
                None) # khu
    

class FV3LMModule(torch.nn.Module):
    def __init__(self, dt, npx, npy, npz, ptop, ak, bk, ntracers, do_phy_mst, do_phy_trb, inputpert_filename):
        super().__init__()

        self.lm = fv3jedi_lm_mod.fv3jedi_lm_type()

        self.lm.conf.inputpert_filename = inputpert_filename

        if do_phy_mst or do_phy_trb:
            self.lm.conf.do_phy = 1
            self.lm.conf.do_phy_trb = 1 if do_phy_trb else 0 
            self.lm.conf.do_phy_mst = 1 if do_phy_mst else 0
        else:
            self.lm.conf.do_phy = 0
            self.lm.conf.do_phy_trb = 0
            self.lm.conf.do_phy_mst = 0

        self.lm.create(dt, npx, npy, npz, ptop, ak, bk)
        
        if ntracers > 0:
            include_pert_tracers = True
            self.lm.allocate_tracers(self.lm.conf.isc, self.lm.conf.iec, self.lm.conf.jsc, self.lm.conf.jec, npz, ntracers, include_pert_tracers)

        self.lm.init_nl()

    def forward(self, u, v, t, delp, tracers, w=None, delz=None,
                cfcn=None, qls=None, qcn=None, phis=None, frocean=None, frland=None, varflt=None, ustar=None,
                bstar=None, zpbl=None, cm=None, ct=None, cq=None, kcbl=None, ts=None, khl=None, khu=None):
        
        return FV3LMFnc.apply(self.lm, u, v, t, delp, tracers, w, delz,
                cfcn, qls, qcn, phis, frocean, frland, varflt, ustar,
                bstar, zpbl, cm, ct, cq, kcbl, ts, khl, khu)

    def extra_repr(self):
        pass

