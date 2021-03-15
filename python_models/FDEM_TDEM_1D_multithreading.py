# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import time
import matplotlib.pyplot as plt
from multiprocessing import Pool

import discretize
from SimPEG import (
    maps,
    utils,
    data_misfit,
    regularization,
    optimization,
    inversion,
    inverse_problem,
    directives,
)
import numpy as np
from SimPEG.electromagnetics import frequency_domain as FDEM, time_domain as TDEM, mu_0
import matplotlib
start = time.perf_counter()

try:
    from pymatsolver import Pardiso as Solver
except ImportError:
    from SimPEG import SolverLU as Solver


def run(sig_layer):

    # Set up cylindrically symmeric mesh
    cs, ncx, ncz, npad = 10.0, 15, 25, 13  # padded cyl mesh
    hx = [(cs, ncx), (cs, npad, 1.3)]
    hz = [(cs, npad, -1.3), (cs, ncz), (cs, npad, 1.3)]
    mesh = discretize.CylMesh([hx, 1, hz], "00C")

    # Conductivity model
    layerz = np.r_[-200.0, -100.0]
    layer = (mesh.vectorCCz >= layerz[0]) & (mesh.vectorCCz <= layerz[1])
    active = mesh.vectorCCz < 0.0
    sig_half = 1e-2  # Half-space conductivity
    sig_air = 1e-8  # Air conductivity
    sigma = np.ones(mesh.nCz) * sig_air
    sigma[active] = sig_half
    sigma[layer] = sig_layer

    # Mapping
    actMap = maps.InjectActiveCells(mesh, active, np.log(1e-8), nC=mesh.nCz)
    mapping = maps.ExpMap(mesh) * maps.SurjectVertical1D(mesh) * actMap
    mtrue = np.log(sigma[active])

    # ----- FDEM problem & survey ----- #
    rxlocs = utils.ndgrid([np.r_[50.0], np.r_[0], np.r_[0.0]])
    bzr = FDEM.Rx.PointMagneticFluxDensitySecondary(rxlocs, "z", "real")
    bzi = FDEM.Rx.PointMagneticFluxDensitySecondary(rxlocs, "z", "imag")

    freqs = np.logspace(2, 3, 5)
    srcLoc = np.array([0.0, 0.0, 0.0])

    print(
        "min skin depth = ",
        500.0 / np.sqrt(freqs.max() * sig_half),
        "max skin depth = ",
        500.0 / np.sqrt(freqs.min() * sig_half),
    )
    print(
        "max x ",
        mesh.vectorCCx.max(),
        "min z ",
        mesh.vectorCCz.min(),
        "max z ",
        mesh.vectorCCz.max(),
    )

    srcList = [
        FDEM.Src.MagDipole([bzr, bzi], freq, srcLoc, orientation="Z") for freq in freqs
    ]

    surveyFD = FDEM.Survey(srcList)
    prbFD = FDEM.Simulation3DMagneticFluxDensity(
        mesh, survey=surveyFD, sigmaMap=mapping, solver=Solver
    )
    rel_err = 0.03
    dataFD = prbFD.make_synthetic_data(mtrue, relative_error=rel_err, add_noise=True)
    dataFD.noise_floor = np.linalg.norm(dataFD.dclean) * 1e-5

    # FDEM inversion
    np.random.seed(1)
    dmisfit = data_misfit.L2DataMisfit(simulation=prbFD, data=dataFD)
    regMesh = discretize.TensorMesh([mesh.hz[mapping.maps[-1].indActive]])
    reg = regularization.Simple(regMesh)
    opt = optimization.InexactGaussNewton(maxIterCG=10)
    invProb = inverse_problem.BaseInvProblem(dmisfit, reg, opt)

    # Inversion Directives
    beta = directives.BetaSchedule(coolingFactor=4, coolingRate=3)
    betaest = directives.BetaEstimate_ByEig(beta0_ratio=2.0)
    target = directives.TargetMisfit()
    directiveList = [beta, betaest, target]

    inv = inversion.BaseInversion(invProb, directiveList=directiveList)
    m0 = np.log(np.ones(mtrue.size) * sig_half)
    reg.alpha_s = 5e-1
    reg.alpha_x = 1.0
    prbFD.counter = opt.counter = utils.Counter()
    opt.remember("xc")
    moptFD = inv.run(m0)

    # TDEM problem
    times = np.logspace(-4, np.log10(2e-3), 10)
    print(
        "min diffusion distance ",
        1.28 * np.sqrt(times.min() / (sig_half * mu_0)),
        "max diffusion distance ",
        1.28 * np.sqrt(times.max() / (sig_half * mu_0)),
    )
    rx = TDEM.Rx.PointMagneticFluxDensity(rxlocs, times, "z")
    src = TDEM.Src.MagDipole(
        [rx],
        waveform=TDEM.Src.StepOffWaveform(),
        loc=srcLoc,  # same src location as FDEM problem
    )

    surveyTD = TDEM.Survey([src])
    prbTD = TDEM.Simulation3DMagneticFluxDensity(
        mesh, survey=surveyTD, sigmaMap=mapping, Solver=Solver
    )
    prbTD.time_steps = [(5e-5, 10), (1e-4, 10), (5e-4, 10)]

    rel_err = 0.03
    dataTD = prbTD.make_synthetic_data(mtrue, relative_error=rel_err, add_noise=True)
    dataTD.noise_floor = np.linalg.norm(dataTD.dclean) * 1e-5

    # TDEM inversion
    dmisfit = data_misfit.L2DataMisfit(simulation=prbTD, data=dataTD)
    regMesh = discretize.TensorMesh([mesh.hz[mapping.maps[-1].indActive]])
    reg = regularization.Simple(regMesh)
    opt = optimization.InexactGaussNewton(maxIterCG=10)
    invProb = inverse_problem.BaseInvProblem(dmisfit, reg, opt)

    # directives
    beta = directives.BetaSchedule(coolingFactor=4, coolingRate=3)
    betaest = directives.BetaEstimate_ByEig(beta0_ratio=2.0)
    target = directives.TargetMisfit()
    directiveList = [beta, betaest, target]

    inv = inversion.BaseInversion(invProb, directiveList=directiveList)
    m0 = np.log(np.ones(mtrue.size) * sig_half)
    reg.alpha_s = 5e-1
    reg.alpha_x = 1.0
    prbTD.counter = opt.counter = utils.Counter()
    opt.remember("xc")
    moptTD = inv.run(m0)

    # Plot the results
    fig, ax = plt.subplots(1, 3, figsize=(10, 12))

    fs = 13  # fontsize
    matplotlib.rcParams["font.size"] = fs

    # Plot the model
    # z_true = np.repeat(mesh.vectorCCz[active][1:], 2, axis=0)
    # z_true = np.r_[mesh.vectorCCz[active][0], z_true, mesh.vectorCCz[active][-1]]
    activeN = mesh.vectorNz <= 0.0 + cs / 2.0
    z_true = np.repeat(mesh.vectorNz[activeN][1:-1], 2, axis=0)
    z_true = np.r_[mesh.vectorNz[activeN][0], z_true, mesh.vectorNz[activeN][-1]]
    sigma_true = np.repeat(sigma[active], 2, axis=0)

    ax[0].semilogx(sigma_true, z_true, "k-", lw=2, label="True")

    ax[0].semilogx(
        np.exp(moptFD),
        mesh.vectorCCz[active],
        "bo",
        ms=6,
        markeredgecolor="k",
        markeredgewidth=0.5,
        label="FDEM",
    )
    ax[0].semilogx(
        np.exp(moptTD),
        mesh.vectorCCz[active],
        "r*",
        ms=10,
        markeredgecolor="k",
        markeredgewidth=0.5,
        label="TDEM",
    )
    ax[0].set_ylim(-700, 0)
    ax[0].set_xlim(5e-3, 1e-1)

    ax[0].set_xlabel("Conductivity (S/m)", fontsize=fs)
    ax[0].set_ylabel("Depth (m)", fontsize=fs)
    ax[0].grid(which="both", color="k", alpha=0.5, linestyle="-", linewidth=0.2)
    ax[0].legend(fontsize=fs, loc=4)

    # plot the data misfits - negative b/c we choose positive to be in the
    # direction of primary

    ax[1].plot(freqs, -dataFD.dobs[::2], "k-", lw=2, label="Obs (real)")
    ax[1].plot(freqs, -dataFD.dobs[1::2], "k--", lw=2, label="Obs (imag)")

    dpredFD = prbFD.dpred(moptTD)
    ax[1].loglog(
        freqs,
        -dpredFD[::2],
        "bo",
        ms=6,
        markeredgecolor="k",
        markeredgewidth=0.5,
        label="Pred (real)",
    )
    ax[1].loglog(
        freqs, -dpredFD[1::2], "b+", ms=10, markeredgewidth=2.0, label="Pred (imag)"
    )

    ax[2].loglog(times, dataTD.dobs, "k-", lw=2, label="Obs")
    ax[2].loglog(
        times,
        prbTD.dpred(moptTD),
        "r*",
        ms=10,
        markeredgecolor="k",
        markeredgewidth=0.5,
        label="Pred",
    )
    ax[2].set_xlim(times.min() - 1e-5, times.max() + 1e-4)

    # Labels, gridlines, etc
    ax[2].grid(which="both", alpha=0.5, linestyle="-", linewidth=0.2)
    ax[1].grid(which="both", alpha=0.5, linestyle="-", linewidth=0.2)

    ax[1].set_xlabel("Frequency (Hz)", fontsize=fs)
    ax[1].set_ylabel("Vertical magnetic field (-T)", fontsize=fs)

    ax[2].set_xlabel("Time (s)", fontsize=fs)
    ax[2].set_ylabel("Vertical magnetic field (T)", fontsize=fs)

    ax[2].legend(fontsize=fs, loc=3)
    ax[1].legend(fontsize=fs, loc=3)
    ax[1].set_xlim(freqs.max() + 1e2, freqs.min() - 1e1)

    ax[0].set_title("(a) Recovered Models", fontsize=fs)
    ax[1].set_title("(b) FDEM observed vs. predicted", fontsize=fs)
    ax[2].set_title("(c) TDEM observed vs. predicted", fontsize=fs)

    plt.tight_layout(pad=1.5)

    return ax

if __name__ == '__main__':
    with Pool(8) as p:
        print(p.map(run, [5e-2, 5e-3, 5e-4]))

finish = time.perf_counter()

print(f'Finished in {round(finish-start, 4)} second(s)')


#y = 2

#def f(x):
#    sleep(2)
#    fig, ax = plt.subplots()
#    ax.scatter(x, y)
#    return ax

