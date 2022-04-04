using Base.Threads
using SparseArrays
using MAT
using DocStringExtensions
using Parameters
using StaticArrays
using BenchmarkTools
using TimerOutputs

const to = TimerOutput()

""""
Static parameters: Grids, markers, switches, constants, etc. which remain
constant throughout the simulation.

$(TYPEDFIELDS)
"""
@with_kw struct StaticParameters
    # Radioactive switches
    "radioactive heating from 26Al active"
    hr_al::Bool = true
    "radioactive heating from 60Fe active"	
    hr_fe::Bool = true
    # Model size and resolution
    "horizontal model size [m]"
    xsize::Int64 = 140000
    "vertical model size [m]"
    ysize::Int64 = 140000
    "horizontal grid resolution"
    Nx::Int
    "vertical grid resolution"	
    Ny::Int
    Nx1::Int64 = Nx + 1
    Ny1::Int64 = Ny + 1
    "horizontal grid step [m]"
    dx::Float64 = xsize / (Nx-1)
    "vertical grid step [m]"
    dy::Float64 = ysize / (Ny-1)
    # Planetary parameters
    "planetary radius [m]"
    rplanet::Int64 = 50000
    "crust radius [m]"
    rcrust::Int64 = 48000
    "surface pressure [Pa]"
    psurface::Float64 = 1e+3
    # Markers
    "number of markers per cell in horizontal direction"
    Nxmc::Int
    "number of markers per cell in vertical direction"
    Nymc::Int
    "marker grid resolution in horizontal direction"
    Nxm::Int64 = (Nx - 1) * Nxmc
    "marker grid resolution in vertical direction"
    Nym::Int64 = (Ny - 1) * Nymc
    "marker grid step in horizontal direction"
    dxm::Float64 = xsize / Nxm
    "marker grid step in vertical direction"
    dym::Float64 = ysize / Nym
    "number of markers at start"
    startmarknum::Int64 = Nxm * Nym
    # Physical constants
    "gravitational constant [m^3*kg^-1*s^-2]"
    G::Float64 = 6.672e-11
    "scaled pressure"    
    pscale::Float64 = 1e+23 / dx
    # Materials properties:         Planet      Crust       Space
    "solid Density [kg/m^3]"
    rhosolidm::SVector{3, Float64}   = [ 3300.0    , 3300.0    ,    1.0    ]
    "fluid density [kg/m^3]"	
    rhofluidm::SVector{3, Float64}   = [ 7000.0    , 7000.0    , 1000.0    ]
    "solid viscosity [Pa*s]"
    etasolidm::SVector{3, Float64}   = [    1.0e+16,    1.0e+16,    1.0e+14]
    "molten solid viscosity [Pa*s]"
    etasolidmm::SVector{3, Float64}  = [    1.0e+14,    1.0e+14,    1.0e+14]
    "fluid viscosity [Pa*s]"
    etafluidm::SVector{3, Float64}   = [    1.0e-02,    1.0e-02,    1.0e+12]
    "molten fluid viscosity [Pa*s]"
    etafluidmm::SVector{3, Float64}  = [    1.0e-02,    1.0e-02,    1.0e+12]
    "solid volumetric heat capacity [kg/m^3]"
    rhocpsolidm::SVector{3, Float64} = [    3.3e+06,    3.3e+06,    3.0e+06]
    "fluid volumetric heat capacity [kg/m^3]"
    rhocpfluidm::SVector{3, Float64} = [    7.0e+06,    7.0e+06,    3.0e+06]
    "solid thermal expansion [1/K]"
    alphasolidm::SVector{3, Float64} = [    3.0e-05,    3.0e-05,    0.0    ]
    "fluid thermal expansion [1/K]"
    alphafluidm::SVector{3, Float64} = [    5.0e-05,    5.0e-05,    0.0    ]
    "solid thermal conductivity [W/m/K]"
    ksolidm::SVector{3, Float64}     = [    3.0    ,    3.0    , 3000.0    ]
    "fluid thermal conductivity [W/m/K]"
    kfluidm::SVector{3, Float64}     = [   50.0    ,   50.0    , 3000.0    ]
    "solid radiogenic heat production [W/m^3]"
    hrsolidm::Array{Float64}    = [    0.0    ,    0.0    ,    0.0    ]
    "fluid radiogenic heat production [W/m^3]"
    hrfluidm::Array{Float64}    = [    0.0    ,    0.0    ,    0.0    ]
    "solid shear modulus [Pa]"
    gggsolidm::SVector{3, Float64}   = [    1.0e+10,    1.0e+10,    1.0e+10]
    "solid friction coefficient"
    frictsolidm::SVector{3, Float64} = [    0.6    ,    0.6    ,    0.0    ]
    "solid compressive strength [Pa]"
    cohessolidm::SVector{3, Float64} = [    1.0e+08,    1.0e+08,    1.0e+08]
    "solid tensile strength [Pa]"
    tenssolidm ::SVector{3, Float64} = [    6.0e+07,    6.0e+07,    6.0e+07]
    "standard permeability [m^2]"
    kphim0::SVector{3, Float64}      = [    1.0e-13,    1.0e-13,    1.0e-17]
    "Coefficient to compute compaction viscosity from shear viscosity"
    etaphikoef::Float64 = 1e-4
    # 26Al decay
    "26Al half life [s]"
    t_half_al::Float64 = 717000 * 31540000
    "26Al decay constant"
    tau_al::Float64 = t_half_al / log(2)
    "initial ratio of 26Al and 27Al isotopes"
    ratio_al::Float64 = 5.0e-5
    "E 26Al [J]"
    E_al::Float64 = 5.0470e-13
    "26Al atoms/kg"
    f_al::Float64 = 1.9e23
    # 60Fe decay
    "60Fe half life [s]"	
    t_half_fe::Float64 = 2620000 * 31540000
    "60Fe decay constant"
    tau_fe::Float64 = t_half_fe / log(2)
    "initial ratio of 60Fe and 56Fe isotopes"	
    ratio_fe::Float64 = 1e-6
    "E 60Fe [J]"	
    E_fe::Float64 = 4.34e-13
    "60Fe atoms/kg"	
    f_fe::Float64 = 1.957e24
    # Melting temperatures
    "silicate melting temperature [K]"
    tmsilicate::Float64 = 1e+6
    "iron melting temperature [K]"
    tmiron::Float64 = 1273 
    # Porosities
    "standard Fe fraction [porosity]"
    phim0::Float64 = 0.2
    "min porosity"	
    phimin::Float64 = 1e-4
    "max porosity"
    phimax::Float64 = 1 - phimin            
    # Mechanical boundary conditions: free slip=-1 / no slip=1
    "mechanical boundary condition left"
    bcleft::Float64 = -1
    "mechanical boundary condition right"
    bcright::Float64 = -1
    "mechanical boundary condition top"
    bctop::Float64 = -1
    "mechanical boundary condition bottom"
    bcbottom::Float64 = -1
    # Hydraulic boundary conditions: free slip=-1 / no slip=1
    "hydraulic boundary condition left"
    bcfleft::Float64 = -1
    "hydraulic boundary condition right"
    bcfright::Float64 = -1
    "hydraulic boundary condition top"
    bcftop::Float64 = -1
    "hydraulic boundary condition bottom"
    bcfbottom::Float64 = -1
    # Extension/shortening velocities
    "shortening strain rate"
    strainrate::Float64 = 0e-13
    "x extension/shortening velocity left"
    vxleft::Float64 = strainrate * xsize / 2
    "x extension/shortening velocity right"
    vxright::Float64= -strainrate * xsize / 2
    "y extension/shortening velocity top"
    vytop::Float64 = - strainrate * ysize / 2
    "y extension/shortening velocity bottom"
    vybottom::Float64 = strainrate * ysize / 2
    # timestepping parameters
    "mat filename"
    nname::String = "madcph_"
    ".mat storage periodicity"
    savematstep::Int64 = 50
    "Maximal computational timestep [s]"
    dtelastic::Float64 = 1e+11 
    "Coefficient to decrease computational timestep"
    dtkoef::Float64 = 2 
    "Coefficient to increase computational timestep"
    dtkoefup::Float64 = 1.1 
    "Number of iterations before changing computational timestep"
    dtstep::Int64 = 200 
    "Max marker movement per time step [grid steps]"
    dxymax::Float64 = 0.05 
    "Weight of averaged velocity for moving markers"
    vpratio::Float64 = 1 / 3 
    "Max temperature change per time step [K]"
    DTmax::Float64 = 20 
    "Subgrid temperature diffusion parameter"
    dsubgridt::Float64 = 0 
    "Subgrid stress diffusion parameter"
    dsubgrids::Float64 = 0
    "length of year [s]"
    yearlength::Float64 = 365.25 * 24 * 3600
    "Time sum (start) [s]"
    starttimesum::Float64 = 1e6 * yearlength 
    "Time sum (end) [s]"
    endtimesum::Float64 = 15 * 1000000 * yearlength
    "Lower viscosity cut-off [Pa s]"	
    etamin::Float64 = 1e+12 
    "Upper viscosity cut-off [Pa s]"
    etamax::Float64 = 1e+23 
    "Number of plastic iterations"
    nplast::Int64 = 100000
    "Periodicity of visualization"
    visstep::Int64 = 1 
    "Tolerance level for yielding error()"
    yerrmax::Float64 = 1e+2 
    "Yielding error of nodes"
    YERRNOD::Array{Float64} = zeros(1, nplast) 
    "Weight for old viscosity"
    etawt::Float64 = 0 
    "max porosity ratio change per time step"
    dphimax::Float64 = 0.01
    "starting timestep"
    startstep::Int64 = 1
    "number of timesteps"
    nsteps::Int64 = 30000 
end


 """
 Dynamic parameters: Grids, markers, switches, constants, etc. which undergo
 change throughout the simulation.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct DynamicParameters
    "timestep counter (current)"
    timestep::Int64
    "computational timestep (current) [s]"
    dt::Float64
    "time sum (current) [s]"
    timesum::Float64
    "current number of markers"
    marknum::Int64
    "inner constructor"
    DynamicParameters(startstep, dt, timesum, marknum) = new(
        startstep, dt, timesum, marknum)
    DynamicParameters(sp::StaticParameters) = new(
        sp.startstep, sp.dt, sp.starttimesum, sp.startmarknum)
    end


"""
Basic node coordinates

$(TYPEDFIELDS)
"""
@with_kw struct BasicNodes
# Base.@kwdef struct BasicNodes
    "horizontal coordinates of basic grid points [m]"
    x::Array{Float64}
    "vertical coordinates of basic grid points [m]"
    y::Array{Float64}
    "inner constructor"
    BasicNodes(xsize, ysize, dx, dy) = new(
        collect(0:dx:xsize),
        collect(0:dy:ysize)
        )
    BasicNodes(sp::StaticParameters) = new(
        collect(0:sp.dx:sp.xsize),
        collect(0:sp.dy:sp.ysize)
        )
end


"""
Vx node coordinates

$(TYPEDFIELDS)
"""
@with_kw struct VxNodes
# Base.@kwdef struct VxNodes
    "horizontal coordinates of vx grid points [m]"
    xvx::Array{Float64}
    "vertical coordinates of vx grid points [m]"
    yvx::Array{Float64}
    "inner constructor"
    VxNodes(xsize, ysize, dx, dy) = new(
        collect(0:dx:xsize+dy),
        collect(-dy/2:dy:ysize+dy/2)
        )
    VxNodes(sp::StaticParameters) = new(
        collect(0:sp.dx:sp.xsize+sp.dy),
        collect(-sp.dy/2:sp.dy:sp.ysize+sp.dy/2)
        )
end


"""
Vy node coordinates

$(TYPEDFIELDS)
"""
@with_kw struct VyNodes
# Base.@kwdef struct VyNodes
    "horizontal coordinates of vy grid points [m]"
    xvy::Array{Float64}
    "vertical coordinates of vy grid points [m]"
    yvy::Array{Float64}
    "inner constructor"
    VyNodes(xsize, ysize, dx, dy) = new(
        collect(-dx/2:dx:xsize+dx/2),
        collect(0:dy:ysize+dy)
        )
    VyNodes(sp::StaticParameters) = new(
        collect(-sp.dx/2:sp.dx:sp.xsize+sp.dx/2),
        collect(0:sp.dy:sp.ysize+sp.dy)
        )
end


"""
P node coordinates

$(TYPEDFIELDS)
"""
@with_kw struct PNodes
# Base.@kwdef struct PNodes
    "horizontal coordinates of P grid points [m]"
    xp::Array{Float64}
    "vertical coordinates of P grid points [m]"
    yp::Array{Float64}
    "inner constructor"
    PNodes(xsize, ysize, dx, dy) = new(
        collect(-dx/2:dx:xsize+dx/2),
        collect(-dy/2:dy:ysize+dy/2)
        )
    PNodes(sp::StaticParameters) = new(
        collect(-sp.dx/2:sp.dx:sp.xsize+sp.dx/2),
        collect(-sp.dy/2:sp.dy:sp.ysize+sp.dy/2)
        )     
end


"""
Basic node properties

$(TYPEDFIELDS)
"""
@with_kw struct BasicNodalArrays
# Base.@kwdef mutable struct BasicNodalArrays
    "viscoplastic viscosity, Pa*s"
    eta::Array{Float64}
    "viscous viscosity, Pa*s"
    eta0::Array{Float64}
    "shear modulus, Pa"
    ggg::Array{Float64}
    "epsilonxy, 1/s"
    exy::Array{Float64}
    "sigma0xy, 1/s"
    sxy0::Array{Float64}
    "rotation rate, 1/s"
    wyx::Array{Float64}
    "compressive strength, Pa"
    coh::Array{Float64}
    "tensile strength, Pa"
    ten::Array{Float64}
    "friction"
    fri::Array{Float64}
    "plastic yielding mark, 1=yes,0=no"
    yny::Array{Int8}
    "inner constructor"
    NodalArrays(Nx, Ny) = new(
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx),
        zeros(Ny,Nx)
        )
    NodalArrays(sp::StaticParameters) = new(
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx),
        zeros(sp.Ny,sp.Nx)
        )
end


"""
Vx node properties

$(TYPEDFIELDS)
"""
@with_kw struct VxNodalArrays
# Base.@kwdef mutable struct VxNodalArrays
    "density [kg/m^3]"
    rhox::Array{Float64}
    "fluid density [kg/m^3]"
    rhofx::Array{Float64}
    "thermal conductivity [W/m/K]"
    kx::Array{Float64}
    "porosity"
    phix::Array{Float64}
    "solid vx-velocity [m/s]"
    vx::Array{Float64}
    "fluid vx-velocity [m/s]"
    vxf::Array{Float64}
    "etafluid/kphi ratio [m^2]"
    rx::Array{Float64}
    "qx-darcy flux [m/s]"
    qxD::Array{Float64}
    "gx-gravity [m/s^2]"
    gx::Array{Float64}
    "inner constructor"
    NodalArrays(Nx1, Ny1) = new(
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1)
        )
    NodalArrays(sp::StaticParameters) = new(
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1)
     )
end


"""
Vy node properties

$(TYPEDFIELDS)
"""
@with_kw struct VyNodalArrays
# Base.@kwdef mutable struct VyNodalArrays
    "density [kg/m^3]"
    rhoy::Array{Float64}
    "fluid density [kg/m^3]"
    rhofy::Array{Float64}
    "thermal conductivity [W/m/K]"
    ky::Array{Float64}
    "porosity"
    phiy::Array{Float64}
    "solid vy-velocity [m/s]"
    vy::Array{Float64}
    "fluid vy-velocity [m/s]"
    vyf::Array{Float64}
    "etafluid/kphi ratio [m^2]"
    ry::Array{Float64}
    "qy-darcy flux [m/s]"
    qyD::Array{Float64}
    "gy-gravity [m/s^2]"
    gy::Array{Float64}
    "inner constructor"
    NodalArrays(Nx1, Ny1) = new(
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1),
        zeros(Ny1,Nx1)
        )
    NodalArrays(sp::StaticParameters) = new(
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1)
     )
end


"""
P node properties

$(TYPEDFIELDS)
"""
@with_kw struct PNodalArrays
# Base.@kwdef mutable struct PNodalArrays
    "density [kg/m^3]"
    rho::Array{Float64}
    "volumetric heat capacity [J/m^3/K]"
    rhocp::Array{Float64}
    "thermal expansion [J/m^3/K]"
    alpha::Array{Float64}
    "fluid thermal expansion [J/m^3/K]"
    alphaf::Array{Float64}
    "radioactive heating [W/m^3]"
    hr::Array{Float64}
    "adiabatic heating [W/m^3]"
    ha::Array{Float64}
    "shear heating [W/m^3]"
    hs::Array{Float64}
    "viscosity [Pa*s]"
    etap::Array{Float64}
    "shear modulus [Pa]"
    gggp::Array{Float64}
    "EPSILONxx [1/s]"
    exx::Array{Float64}
    "SIGMA'xx [1/s]"
    sxx::Array{Float64}
    "SIGMA0'xx [1/s]"
    sxx0::Array{Float64}
    "old temperature [K]"
    tk1::Array{Float64}
    "new temperature [K]"
    tk2::Array{Float64}
    "solid Vx in pressure nodes [m/s]"
    vxp::Array{Float64}
    "solid Vy in pressure nodes [m/s]"
    vyp::Array{Float64}
    "fluid Vx in pressure nodes [m/s]"
    vxpf::Array{Float64}
    "fluid Vy in pressure nodes [m/s]"
    vypf::Array{Float64}
    "total pressure [Pa]"
    pr::Array{Float64}
    "fluid pressure [Pa]"
    pf::Array{Float64}
    "solid pressure [Pa]"
    ps::Array{Float64}
    "old total pressure [Pa]"
    pr0::Array{Float64}
    "old fluid pressure [Pa]"
    pf0::Array{Float64}
    "old solid pressure [Pa]"
    ps0::Array{Float64}
    "bulk viscosity [Pa*s]"
    etaphi::Array{Float64}
    "bulk compresibility [Pa*s]"
    bettaphi::Array{Float64}
    "porosity"
    phi::Array{Float64}
    "Dln[(1-PHI)/PHI]/Dt"
    aphi::Array{Float64}
    "gravity potential [J/kg]"
    fi::Array{Float64}
    "inner constructor"
    PNodalArrays(Nx1, Ny1) = new(
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1),
        zeros(Ny1, Nx1)
    )
    PNodalArrays(sp::StaticParameters) = new(
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1),
        zeros(sp.Ny1, sp.Nx1)
    )
end


"""
Marker properties: Fixed and calculated during timestepping

$(TYPEDFIELDS)
"""
@with_kw struct MarkerArrays
    # original marker properties 
    "horizontal coordinates [m]"
    xm::Vector{Float64}
    "vertical coordinates [m]"
    ym::Vector{Float64}
    "material type"
    tm::Vector{Int8}
    "marker temperature [K]"
    tkm::Vector{Float64}
    "SIGMA'xx [Pa]"
    sxxm::Vector{Float64}
    "SIGMAxy [Pa]"
    sxym::Vector{Float64}
    "Visco-plastic viscosity [Pa]"
    etavpm::Vector{Float64}
    "Marker porosity"
    phim::Vector{Float64}
    # fixed marker properties used during timestepping calculations
    # N/A: omitted - reconsider?
    # marker properties calculated during timestepping
    "kphim"
    kphim::Vector{Float64}
    "rhototalm"
    rhototalm::Vector{Float64}
    "rhocptotalm"
    rhocptotalm::Vector{Float64}
    "etatotalm"
    etatotalm::Vector{Float64}
    "hrtotalm"
    hrtotalm::Vector{Float64}
    "ktotalm"
    ktotalm::Vector{Float64}
    "gggtotalm"
    gggtotalm::Vector{Float64}
    "fricttotalm"
    fricttotalm::Vector{Float64}
    "cohestotalm"
    cohestotalm::Vector{Float64}
    "tenstotalm"
    tenstotalm::Vector{Float64}
    "etafluidcur"
    etafluidcur::Vector{Float64}
    "rhofluidcur"
    rhofluidcur::Vector{Float64}
    "inner constructor"
    MarkerArrays(marknum) = new(
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum)
    )
end


"""
Interpolation arrays: Calculated during timestepping
to interpolate properties from markers to nodes

$(TYPEDFIELDS)
"""
@with_kw struct InterpArrays
    # basic nodes
    "basic nodes: ETA0SUM"
    ETA0SUM::Tuple
    "basic nodes: ETASUM"
    ETASUM::Tuple
    "basic nodes: GGGSUM"
    GGGSUM::Tuple
    "basic nodes: SXYSUM"
    SXYSUM::Tuple
    "basic nodes: COHSUM"
    COHSUM::Tuple
    "basic nodes: TENSUM"
    TENSUM::Tuple
    "basic nodes: FRISUM"	
    FRISUM::Tuple
    "basic nodes: WTSUM"
    WTSUM::Tuple
    # Vx-nodes
    "Vx-nodes: RHOXSUM"	
    RHOXSUM::Tuple
    "Vx-nodes: RHOFXSUM"
    RHOFXSUM::Tuple
    "Vx-nodes: KXSUM"
    KXSUM::Tuple
    "Vx-nodes: PHIXSUM"
    PHIXSUM::Tuple
    "Vx-nodes: RXSUM"	
    RXSUM::Tuple
    "Vx-nodes: WTXSUM"	
    WTXSUM::Tuple
    # Vy-nodes
    "Vy-nodes: RHOYSUM"	
    RHOYSUM::Tuple
    "Vy-nodes: RHOFYSUM"
    RHOFYSUM::Tuple
    "Vy-nodes: KYSUM"	
    KYSUM::Tuple
    "Vy-nodes: PHIYSUM"
    PHIYSUM::Tuple
    "Vy-nodes: RYSUM"
    RYSUM::Tuple
    "Vy-nodes: WTYSUM"
    WTYSUM::Tuple
    # P-Nodes
    "P-nodes: GGGPSUM"
    GGGPSUM::Tuple
    "P-nodes: SXXSUM"
    SXXSUM::Tuple
    "P-nodes: RHOSUM"
    RHOSUM::Tuple
    "P-nodes: RHOCPSUM"
    RHOCPSUM::Tuple
    "P-nodes: ALPHASUM"	
    ALPHASUM::Tuple
    "P-nodes: ALFAFSUM"
    ALPHAFSUM::Tuple
    "P-nodes: HRSUM"	
    HRSUM::Tuple
    "P-nodes: TKSUM"
    TKSUM::Tuple
    "P-nodes: PHISUM"
    PHISUM::Tuple
    "P-nodes: WTPSUM"
    WTPSUM::Tuple
    "inner constructor"
    InterpArrays(Nx, Ny, Nx1, Ny1) = new(
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny, Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, Ny1, Nx1) for _ in 1:nthreads()]),
    )
    InterpArrays(sp::StaticParameters) = new(
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny, sp.Nx) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()]),
        Tuple([Matrix{Float64}(undef, sp.Ny1, sp.Nx1) for _ in 1:nthreads()])
    )
end


"""
Global matrices: Hydro-mechanical solution

$(TYPEDFIELDS)
"""
@with_kw struct GlobalHydroMechanicalSolution
# Base.@kwdef mutable struct GlobalHydroMechanicalSolution
    "L matrix"
    L::SparseMatrixCSC{Float64, Int64}
    "R vector"
    R::Vector{Float64}
    "inner constructor"
    GlobalHydroMechanical(Nx1, Ny1) = new(
        spzeros(Nx1*Ny1*6, Nx1*Ny1*6),
        zeros(Nx1*Ny1*6)
        )
    GlobalHydroMechanical(sp::StaticParameters) = new(
        spzeros(sp.Nx1*sp.Ny1*6, sp.Nx1*sp.Ny1*6),
        zeros(sp.Nx1*sp.Ny1*6)
        )
end


"""
Global matrices: Thermal solution

$(TYPEDFIELDS)
"""
@with_kw struct GlobalThermalSolution
# Base.@kwdef mutable struct GlobalThermalSolution
    "LT matrix"
    LT::SparseMatrixCSC{Float64, Int64}
    "RT vector"
    RT::Vector{Float64}
    "inner constructor"
    GlobalThermal(Nx1, Ny1) = new(
        spzeros(Nx1*Ny1, Nx1*Ny1),
        zeros(Nx1*Ny1)
        )
    GlobalThermal(sp::StaticParameters) = new(
        spzeros(sp.Nx1*sp.Ny1, sp.Nx1*sp.Ny1),
        zeros(sp.Nx1*sp.Ny1)
        )
end


"""
Global matrices: Gravity solution

$(TYPEDFIELDS)
"""
@with_kw struct GlobalGravitySolution
# Base.@kwdef mutable struct GlobalGravitySolution
    "LP matrix"
    LP::SparseMatrixCSC{Float64, Int64}
    "RP vector"
    RP::Vector{Float64}
    "inner constructor"
    GlobalGravity(Nx1, Ny1) = new(
        spzeros(Nx1*Ny1, Nx1*Ny1),
        zeros(Nx1*Ny1)
        )
    GlobalGravity(sp::StaticParameters) = new(
        spzeros(sp.Nx1*sp.Ny1, sp.Nx1*sp.Ny1),
        zeros(sp.Nx1*sp.Ny1)
        )
end



"""
Initialize markers according to model parameters

$(SIGNATURES)

# Details

    - ma: arrays containing marker properties
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

# Returns

    - nothing
"""
function initmarkers!(
    ma::MarkerArrays, sp::StaticParameters, dp::DynamicParameters)
    @unpack xm, ym, tm, tkm, phim, etavpm = ma
    @unpack xsize, ysize, Nxm, Nym, dxm, dym, rplanet, rcrust, phimin,
        etasolidm, phim0 = sp

    xcenter = xsize / 2
    ycenter = ysize / 2
    radius(x, y) = sqrt((x - xcenter)^2 + (y - ycenter)^2)
    for jm=1:1:Nxm, im=1:1:Nym
        # calculate marker counter
        m = (jm-1) * Nym + im
        # define marker coordinates
        xm[m] = dxm/2 + (jm-1) * dxm + (rand()-0.5) * dxm
        ym[m] = dym/2 + (im-1) * dym + (rand()-0.5) * dym
        # primary marker properties 
        rmark = radius(xm[m], ym[m])
        if rmark < rplanet
            # planet
            if rmark > rcrust
                # crust
                tm[m] = 2
            else
                # mantle
                tm[m] = 1
            end
            # temperature
            tkm[m] = 300
            # porosity
            phim[m] = phim0 * (1.0 + (rand()-0.5))
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]] # *exp(-28*phim[m])
        else
            # sticky space ("air") [to have internal free surface]
            # space
            tm[m] = 3
            # temperature
            tkm[m] = 273
            # porosity
            phim[m] = phimin
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]]
        end
        # secondary marker properties
        compute_static_marker_params!(m, ma, sp)
        compute_dynamic_marker_params!(m, ma, sp, dp)
    end
end


"""
Compute static marker properties.
Runs once during initialization.

$(SIGNATURES)

# Details 

    - m: marker counter of marker whose static properties are to be computed
    - ma: arrays containing marker properties
    - sp: static simulation parameters

# Returns

    - nothing
"""	
function compute_static_marker_params!(
    m::Int64, ma::MarkerArrays, sp::StaticParameters)
    # @unpack_MarkerArrays ma
    @unpack tm, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, gggtotalm,
        fricttotalm, cohestotalm, tenstotalm, etafluidcur, rhofluidcur = ma
    @unpack rhosolidm, rhocpsolidm, etasolidm, hrsolidm, ksolidm, gggsolidm,
        frictsolidm, cohessolidm, tenssolidm, etafluidm, rhofluidm = sp

    # static secondary marker properties
    if tm[m] < 3
        # rocks
        # pass
    else
        # air
        rhototalm[m] = rhosolidm[tm[m]]
        rhocptotalm[m] = rhocpsolidm[tm[m]]
        etatotalm[m] = etasolidm[tm[m]]
        hrtotalm[m] = hrsolidm[tm[m]]
        ktotalm[m] = ksolidm[tm[m]]
    end
    # common for rocks and air
    gggtotalm[m] = gggsolidm[tm[m]]
    fricttotalm[m] = frictsolidm[tm[m]]
    cohestotalm[m] = cohessolidm[tm[m]]
    tenstotalm[m] = tenssolidm[tm[m]]
    etafluidcur[m] = etafluidm[tm[m]]
    rhofluidcur[m] = rhofluidm[tm[m]]
end


"""
Compute convex combination of fluid and solid properties to get total property.

$(SIGNATURES)

# Details

    - fluid: fluid properties
    - solid: solid properties
    - phi: fraction of fluid

# Returns

    - total: computed total property
"""
function total(solid, fluid, phi)
    return solid * (1.0-phi) + fluid * phi
end


"""
Compute temperature-dependent total viscosity for iron-containing silicate rock.

$(SIGNATURES)

# Details

    - tk: temperature [K]
    - tmsilicate: melting temperature of silicate in [K]
    - tmiron: melting temperature of iron in K [K]
    - etamin: minimum viscosity [Pa s]
    - etasolidm: solid viscosity [Pa s]
    - etasolidmm: molten solid viscosity [Pa s]
    - etafluidm: fluid viscosity [Pa s]
    - etafluidmm: molten fluid viscosity [Pa s]

# Returns

    - etatotal: temperature-dependent total viscosity [Pa s]
"""	
function etatotal_rock(
    tk,
    tmsilicate,
    tmiron,
    etamin,
    etasolidm,
    etasolidmm,
    etafluidm,
    etafluidmm
    )
    return max(
        etamin,
        tk > tmsilicate ? etasolidmm : etasolidm,
        tk > tmiron ? etafluidmm : etafluidm
        )
end

"""
Compute total thermal conductivity of two-phase material.

$(SIGNATURES)

# Details

    - ksolid: solid thermal conductivity [W/m/K]
    - kfluid: fluid thermal conductivity [W/m/K]
    - phi: fraction of solid

# Returns

    - ktotal: total thermal conductivity of mixed phase [W/m/K]
"""
function ktotal(ksolid, kfluid, phi)
    return (ksolid * kfluid/2 + ((ksolid * (3*phi-2)
                                 + kfluid * (1.0-3.0*phi))^2)/16)^0.5
            - (ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))/4
end


"""
Compute iron porosity-dependent permeability.

$(SIGNATURES)

# Details

    - kphim0: standard permeability [m^2]
    - phim: actual (marker) porosity
    - phim0: standard iron fraction (porosity)

# Returns

    - kphim: iron porosity-dependent permeability [m^2]
"""
function kphi(kphim0, phim, phim0)
    return kphim0 * (phim/phim0)^3 / ((1.0-phim)/(1.0-phim0))^2
end


"""
Compute dynamic marker properties.
Runs every time step.

$(SIGNATURES)

# Detail

    - m: marker index of marker whose parameters are to be computed
    - ma: marker arrays containing marker properties to be updated
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

# Returns
    - nothing
"""
function compute_dynamic_marker_params!(
    m::Int64, ma::MarkerArrays, sp::StaticParameters, dp::DynamicParameters
    )
    # @unpack_MarkerArrays ma
    @unpack tm, tkm, phim, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm,
         kphim = ma
    # @unpack_Params sp
    @unpack rhosolidm, rhofluidm, rhocpsolidm, rhocpfluidm, tmiron, tmsilicate,
        etamin, etasolidmm, etasolidm, etafluidmm, etafluidm, kphim0, phim0,
        hrsolidm, hrfluidm, ksolidm, kfluidm  = sp

    if tm[m] < 3
        # rocks
        rhototalm[m] = total(rhosolidm[tm[m]], rhofluidm[tm[m]], phim[m])
        rhocptotalm[m] = total(rhocpsolidm[tm[m]], rhocpfluidm[tm[m]], phim[m])
        etatotalm[m] = etatotal_rock(
            tkm[m],
            tmsilicate,
            tmiron,
            etamin,
            etasolidm[tm[m]],
            etasolidmm[tm[m]],
            etafluidm[tm[m]],
            etafluidmm[tm[m]]
            )
        hrtotalm[m] = total(hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
        ktotalm[m] = ktotal(ksolidm[tm[m]], kfluidm[tm[m]], phim[m])
        
    else
        # air
    
    end
    # common for rocks and air
    kphim[m] = kphi(kphim0[tm[m]], phim0, phim[m])
end


"""
Compute radiogenic heat production of isotope mixture.

$(SIGNATURES)

# Details

    - f: fraction of radioactive matter [atoms/kg]
    - ratio: initial ratio of radioactive to non-radioactive isotopes
    - E: heat energy [J]
    - tau: exponential decay mean lifetime `\tau = \frac{t_{1/2}}{\log{2}}` [s]
    - time: time elapsed since start of radioactive decay [s]

# Returns

    - Q: radiogenic heat production [W/kg]
"""
function Q_radiogenic(f, ratio, E, tau, time)
    return f * ratio * E * exp(-time/tau) / tau
end


"""
Compute radiogenic heat production of 26Al and 60Fe isotopes.

$(SIGNATURES)

# Details
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

# Returns
    - hrsolidm: radiogenic heat production of 26Al [W/m^3]
    - hrfluidm: radiogenic heat production of 60Fe [W/m^3]
"""
function calculate_radioactive_heating(
    sp::StaticParameters, dp::DynamicParameters
    )
    @unpack hr_al, f_al, ratio_al, E_al, tau_al, hr_fe, f_fe, ratio_fe, E_fe,
        tau_fe, rhosolidm, rhofluidm = sp
    @unpack timesum = dp
    #26Al: planet ✓, crust ✓, space ×
    if hr_al
        # 26Al radiogenic heat production [W/kg]
        Q_al = Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum)
        # Solid phase 26Al radiogenic heat production [W/m^3]
        hrsolidm = @SVector [Q_al*rhosolidm[1], Q_al*rhosolidm[2], 0.0]
    else
        hrsolidm = @SVector [0.0, 0.0, 0.0]
    end    
    #60Fe: planet ✓, crust ×, space ×
    if hr_fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Fluid phase 60Fe radiogenic heat production [W/m^3]
        hrfluidm = @SVector [Q_fe*rhofluidm[1], 0.0, 0.0]
    else
        hrfluidm = @SVector [0.0, 0.0, 0.0]
    end
    return hrsolidm, hrfluidm
end


"""
Reset interpolation arrays for next time step.
Avoids costly realllocation of arrays.

$(SIGNATURES)

# Detail

    - interp_arrays: interpolation arrays for current time step

# Returns

    - nothing
"""
function reset_interpolation_arrays!(interp_arrays::InterpArrays)
    @unpack_InterpArrays interp_arrays
    for threadid = 1:1:nthreads()
        # basic nodes
        ETA0SUM[threadid] .= zero(0.0)
        ETASUM[threadid] .= zero(0.0)
        GGGSUM[threadid] .= zero(0.0)
        SXYSUM[threadid] .= zero(0.0)
        COHSUM[threadid] .= zero(0.0)
        TENSUM[threadid] .= zero(0.0)
        FRISUM[threadid] .= zero(0.0)
        WTSUM[threadid] .= zero(0.0)
        # Vx-nodes
        RHOXSUM[threadid] .= zero(0.0)
        RHOFXSUM[threadid] .= zero(0.0)
        KXSUM[threadid] .= zero(0.0)
        PHIXSUM[threadid] .= zero(0.0)
        RXSUM[threadid] .= zero(0.0)
        WTXSUM[threadid] .= zero(0.0)
        # Vy-nodes
        RHOYSUM[threadid] .= zero(0.0)
        RHOFYSUM[threadid] .= zero(0.0)
        KYSUM[threadid] .= zero(0.0)
        PHIYSUM[threadid] .= zero(0.0)
        RYSUM[threadid] .= zero(0.0)
        WTYSUM[threadid] .= zero(0.0)
        # P-Nodes
        GGGPSUM[threadid] .= zero(0.0)
        SXXSUM[threadid] .= zero(0.0)
        RHOSUM[threadid] .= zero(0.0)
        RHOCPSUM[threadid] .= zero(0.0)
        ALPHASUM[threadid] .= zero(0.0)
        ALPHAFSUM[threadid] .= zero(0.0)
        HRSUM[threadid] .= zero(0.0)
        TKSUM[threadid] .= zero(0.0)
        PHISUM[threadid] .= zero(0.0)
        WTPSUM[threadid] .= zero(0.0)
    end
end


# initialize parameters
const static_parameters = StaticParameters(
    Nx = 141,
    Ny = 141,
    Nxmc = 4,
    Nymc = 4
)
# initialize dynamic parameters from static parameters
const dynamic_parameters = DynamicParameters(static_parameters)
const basicnodes = BasicNodes(para.xsize, para.ysize, para.dx, para.dy)
const vxnodes = VxNodes(para.xsize, para.ysize, para.dx, para.dy)
const vynodes = VyNodes(para.xsize, para.ysize, para.dx, para.dy)
const pnodes = PNodes(para.xsize, para.ysize, para.dx, para.dy)

const markers = MarkerArrays(para.marknum)
initmarkers!(markers, static_parameters, dynamic_parameters)

"""
Main simulation loop: runs timestepping.

$(SIGNATURES)

# Detail

    - markers: arrays containing all marker properties
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

# Returns
    
    - nothing
"""
function timestepping(
    markers::MarkerArrays, sp::StaticParameters, dp::DynamicParameters
    )
    # unpack simulation parameters
    @unpack Nx, Ny, Nx1, Ny1, startstep, dtelastic, starttimesum = sp
    @unpack timestep, dtelastic, timesum, marknum, hrsolidm, hrfluidm = dp

    # initialize counters and timestepping loop variables
    "timestep counter (current)"
    timestep::Int64 = startstep
    "computational timestep (current) [s]"
    dt::Float64 = dtelastic
    "time sum (current) [s]"
    timesum::Float64 = starttimesum

    # set up marker interpolation arrays
    interp_arrays = InterpArrays(Nx, Ny, Nx1, Ny1)

    # iterate timesteps   
    for timestep = timestep:1:1000#nsteps
        # set interpolation arrays to zero for this timestep
        reset_interpolation_arrays!(interp_arrays)        

        # calculate radioactive heating
        hrsolidm, hrfluidm = calculate_radioactive_heating(sp, dp)

        # compute marker parameters 
        @threads for m=1:1:marknum
            compute_dynamic_marker_params!(m, markers, sp, dp)
        end

        if timestep % 100 == 0
            println("timestep: ", timestep)
        end

        if timesum > endtimesum
            break
        end

    end # for timestep = timestep:1:nsteps
end # function timestepping(p::Params)



# Save old stresses - RMK: not used anywhere in code
# sxxm00=sxxm; 
# sxym00=sxym;    
    
# # Interpolate properties from markers to nodes
# # Basic nodes
# ETA0SUM=zeros(Ny,Nx)
# ETASUM=zeros(Ny,Nx)
# GGGSUM=zeros(Ny,Nx)
# SXYSUM=zeros(Ny,Nx)
# COHSUM=zeros(Ny,Nx)
# TENSUM=zeros(Ny,Nx)
# FRISUM=zeros(Ny,Nx)
# WTSUM=zeros(Ny,Nx)
# # Vx-nodes
# RHOXSUM=zeros(Ny1,Nx1)
# RHOFXSUM=zeros(Ny1,Nx1)
# KXSUM=zeros(Ny1,Nx1)
# PHIXSUM=zeros(Ny1,Nx1)
# RXSUM=zeros(Ny1,Nx1)
# WTXSUM=zeros(Ny1,Nx1)
# # Vy-nodes
# RHOYSUM=zeros(Ny1,Nx1)
# RHOFYSUM=zeros(Ny1,Nx1)
# KYSUM=zeros(Ny1,Nx1)
# PHIYSUM=zeros(Ny1,Nx1)
# RYSUM=zeros(Ny1,Nx1)
# WTYSUM=zeros(Ny1,Nx1)
# # P-Nodes
# GGGPSUM=zeros(Ny1,Nx1)
# SXXSUM=zeros(Ny1,Nx1)
# RHOSUM=zeros(Ny1,Nx1)
# RHOCPSUM=zeros(Ny1,Nx1)
# ALPHASUM=zeros(Ny1,Nx1)
# ALPHAFSUM=zeros(Ny1,Nx1)
# HRSUM=zeros(Ny1,Nx1)
# TKSUM=zeros(Ny1,Nx1)
# PHISUM=zeros(Ny1,Nx1)
# WTPSUM=zeros(Ny1,Nx1)

for m=1:1:marknum
        # Compute marker parameters
    if(tm[m]<3)
        # Rocks
        kphim=kphim0[tm[m]]*(phim[m]/phim0)^3/((1-phim[m])/(1-phim0))^2; #Permeability
        rhototalm=rhosolidm[tm[m]]*(1-phim[m])+rhofluidm[tm[m]]*phim[m]
        rhocptotalm=rhocpsolidm[tm[m]]*(1-phim[m])+rhocpfluidm[tm[m]]*phim[m]
        etasolidcur=etasolidm[tm[m]]
        if(tkm[m]>tmsilicate)
            etasolidcur=etasolidmm[tm[m]]
        end
        hrtotalm=hrsolidm[tm[m]]*(1-phim[m])+hrfluidm[tm[m]]*phim[m]
        ktotalm=(ksolidm[tm[m]]*kfluidm[tm[m]]/2+((ksolidm[tm[m]]*(3*phim[m]-2)+kfluidm[tm[m]]*(1-3*phim[m]))^2)/16)^0.5-(ksolidm[tm[m]]*(3*phim[m]-2)+ kfluidm[tm[m]]*(1-3*phim[m]))/4
        gggtotalm=gggsolidm[tm[m]]
        fricttotalm=frictsolidm[tm[m]]
        cohestotalm=cohessolidm[tm[m]]
        tenstotalm=tenssolidm[tm[m]]
        etafluidcur=etafluidm[tm[m]]
        rhofluidcur=rhofluidm[tm[m]]
        if(tkm[m]>tmiron)
            etafluidcur=etafluidmm[tm[m]]
        end
        etatotalm=max(etamin,maximum(etafluidcur,etasolidcur));#*exp(-28*phim[m])))
    else()
        # Sticky air
        kphim=kphim0[tm[m]]*(phim[m]/phim0)^3/((1-phim[m])/(1-phim0))^2; #Permeability
        rhototalm=rhosolidm[tm[m]]
        rhocptotalm=rhocpsolidm[tm[m]]
        etatotalm=etasolidm[tm[m]]
        hrtotalm=hrsolidm[tm[m]]
        ktotalm=ksolidm[tm[m]]
        gggtotalm=gggsolidm[tm[m]]
        fricttotalm=frictsolidm[tm[m]]
        cohestotalm=cohessolidm[tm[m]]
        tenstotalm=tenssolidm[tm[m]]
        rhofluidcur=rhofluidm[tm[m]]
        etafluidcur=etafluidm[tm[m]]
    end

  
    # Interpolation to basic nodes
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-x[1])/dx)+1
    i=trunc((ym[m]-y[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<1)
        i=1
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-x[j]
    dymi=ym[m]-y[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update properties
    # i;j Node
    ETA0SUM[i,j]=ETA0SUM[i,j]+etatotalm*wtmij
    ETASUM[i,j]=ETASUM[i,j]+etavpm[m]*wtmij
    GGGSUM[i,j]=GGGSUM[i,j]+1/gggtotalm*wtmij
    SXYSUM[i,j]=SXYSUM[i,j]+sxym[m]*wtmij
    COHSUM[i,j]=COHSUM[i,j]+cohestotalm*wtmij
    TENSUM[i,j]=TENSUM[i,j]+tenstotalm*wtmij
    FRISUM[i,j]=FRISUM[i,j]+fricttotalm*wtmij
    WTSUM[i,j]=WTSUM[i,j]+wtmij
    # i+1;j Node
    ETA0SUM[i+1,j]=ETA0SUM[i+1,j]+etatotalm*wtmi1j
    ETASUM[i+1,j]=ETASUM[i+1,j]+etavpm[m]*wtmi1j
    GGGSUM[i+1,j]=GGGSUM[i+1,j]+1/gggtotalm*wtmi1j
    SXYSUM[i+1,j]=SXYSUM[i+1,j]+sxym[m]*wtmi1j
    COHSUM[i+1,j]=COHSUM[i+1,j]+cohestotalm*wtmi1j
    TENSUM[i+1,j]=TENSUM[i+1,j]+tenstotalm*wtmi1j
    FRISUM[i+1,j]=FRISUM[i+1,j]+fricttotalm*wtmi1j
    WTSUM[i+1,j]=WTSUM[i+1,j]+wtmi1j
    # i;j+1 Node
    ETA0SUM[i,j+1]=ETA0SUM[i,j+1]+etatotalm*wtmij1
    ETASUM[i,j+1]=ETASUM[i,j+1]+etavpm[m]*wtmij1
    GGGSUM[i,j+1]=GGGSUM[i,j+1]+1/gggtotalm*wtmij1
    SXYSUM[i,j+1]=SXYSUM[i,j+1]+sxym[m]*wtmij1
    COHSUM[i,j+1]=COHSUM[i,j+1]+cohestotalm*wtmij1
    TENSUM[i,j+1]=TENSUM[i,j+1]+tenstotalm*wtmij1
    FRISUM[i,j+1]=FRISUM[i,j+1]+fricttotalm*wtmij1
    WTSUM[i,j+1]=WTSUM[i,j+1]+wtmij1
    # i+1;j+1 Node
    ETA0SUM[i+1,j+1]=ETA0SUM[i+1,j+1]+etatotalm*wtmi1j1
    ETASUM[i+1,j+1]=ETASUM[i+1,j+1]+etavpm[m]*wtmi1j1
    GGGSUM[i+1,j+1]=GGGSUM[i+1,j+1]+1/gggtotalm*wtmi1j1
    SXYSUM[i+1,j+1]=SXYSUM[i+1,j+1]+sxym[m]*wtmi1j1
    COHSUM[i+1,j+1]=COHSUM[i+1,j+1]+cohestotalm*wtmi1j1
    TENSUM[i+1,j+1]=TENSUM[i+1,j+1]+tenstotalm*wtmi1j1
    FRISUM[i+1,j+1]=FRISUM[i+1,j+1]+fricttotalm*wtmi1j1
    WTSUM[i+1,j+1]=WTSUM[i+1,j+1]+wtmi1j1;    
    
    
    # Interpolation to vx-nodes
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xvx[1])/dx)+1
    i=trunc((ym[m]-yvx[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xm[m]-xvx[j]
    dymi=ym[m]-yvx[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update properties
    # i;j Node
    RHOXSUM[i,j]=RHOXSUM[i,j]+rhototalm*wtmij
    RHOFXSUM[i,j]=RHOFXSUM[i,j]+rhofluidcur*wtmij
    KXSUM[i,j]=KXSUM[i,j]+ktotalm*wtmij
    PHIXSUM[i,j]=PHIXSUM[i,j]+phim[m]*wtmij
    RXSUM[i,j]=RXSUM[i,j]+etafluidcur/kphim*wtmij
    WTXSUM[i,j]=WTXSUM[i,j]+wtmij
    # i+1;j Node
    RHOXSUM[i+1,j]=RHOXSUM[i+1,j]+rhototalm*wtmi1j
    RHOFXSUM[i+1,j]=RHOFXSUM[i+1,j]+rhofluidcur*wtmi1j
    KXSUM[i+1,j]=KXSUM[i+1,j]+ktotalm*wtmi1j
    PHIXSUM[i+1,j]=PHIXSUM[i+1,j]+phim[m]*wtmi1j
    RXSUM[i+1,j]=RXSUM[i+1,j]+etafluidcur/kphim*wtmi1j
    WTXSUM[i+1,j]=WTXSUM[i+1,j]+wtmi1j
    # i;j+1 Node
    RHOXSUM[i,j+1]=RHOXSUM[i,j+1]+rhototalm*wtmij1
    RHOFXSUM[i,j+1]=RHOFXSUM[i,j+1]+rhofluidcur*wtmij1
    KXSUM[i,j+1]=KXSUM[i,j+1]+ktotalm*wtmij1
    PHIXSUM[i,j+1]=PHIXSUM[i,j+1]+phim[m]*wtmij1
    RXSUM[i,j+1]=RXSUM[i,j+1]+etafluidcur/kphim*wtmij1
    WTXSUM[i,j+1]=WTXSUM[i,j+1]+wtmij1
    # i+1;j+1 Node
    RHOXSUM[i+1,j+1]=RHOXSUM[i+1,j+1]+rhototalm*wtmi1j1
    RHOFXSUM[i+1,j+1]=RHOFXSUM[i+1,j+1]+rhofluidcur*wtmi1j1
    KXSUM[i+1,j+1]=KXSUM[i+1,j+1]+ktotalm*wtmi1j1
    PHIXSUM[i+1,j+1]=PHIXSUM[i+1,j+1]+phim[m]*wtmi1j1
    RXSUM[i+1,j+1]=RXSUM[i+1,j+1]+etafluidcur/kphim*wtmi1j1
    WTXSUM[i+1,j+1]=WTXSUM[i+1,j+1]+wtmi1j1

    # Interpolation to vy-nodes
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xvy[1])/dx)+1
    i=trunc((ym[m]-yvy[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-xvy[j]
    dymi=ym[m]-yvy[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update properties
    # i;j Node
    RHOYSUM[i,j]=RHOYSUM[i,j]+rhototalm*wtmij
    RHOFYSUM[i,j]=RHOFYSUM[i,j]+rhofluidcur*wtmij
    KYSUM[i,j]=KYSUM[i,j]+ktotalm*wtmij
    PHIYSUM[i,j]=PHIYSUM[i,j]+phim[m]*wtmij
    RYSUM[i,j]=RYSUM[i,j]+etafluidcur/kphim*wtmij
    WTYSUM[i,j]=WTYSUM[i,j]+wtmij
    # i+1;j Node
    RHOYSUM[i+1,j]=RHOYSUM[i+1,j]+rhototalm*wtmi1j
    RHOFYSUM[i+1,j]=RHOFYSUM[i+1,j]+rhofluidcur*wtmi1j
    KYSUM[i+1,j]=KYSUM[i+1,j]+ktotalm*wtmi1j
    PHIYSUM[i+1,j]=PHIYSUM[i+1,j]+phim[m]*wtmi1j
    RYSUM[i+1,j]=RYSUM[i+1,j]+etafluidcur/kphim*wtmi1j
    WTYSUM[i+1,j]=WTYSUM[i+1,j]+wtmi1j
    # i;j+1 Node
    RHOYSUM[i,j+1]=RHOYSUM[i,j+1]+rhototalm*wtmij1
    RHOFYSUM[i,j+1]=RHOFYSUM[i,j+1]+rhofluidcur*wtmij1
    KYSUM[i,j+1]=KYSUM[i,j+1]+ktotalm*wtmij1
    PHIYSUM[i,j+1]=PHIYSUM[i,j+1]+phim[m]*wtmij1
    RYSUM[i,j+1]=RYSUM[i,j+1]+etafluidcur/kphim*wtmij1
    WTYSUM[i,j+1]=WTYSUM[i,j+1]+wtmij1
    # i+1;j+1 Node
    RHOYSUM[i+1,j+1]=RHOYSUM[i+1,j+1]+rhototalm*wtmi1j1
    RHOFYSUM[i+1,j+1]=RHOFYSUM[i+1,j+1]+rhofluidcur*wtmi1j1
    KYSUM[i+1,j+1]=KYSUM[i+1,j+1]+ktotalm*wtmi1j1
    PHIYSUM[i+1,j+1]=PHIYSUM[i+1,j+1]+phim[m]*wtmi1j1
    RYSUM[i+1,j+1]=RYSUM[i+1,j+1]+etafluidcur/kphim*wtmi1j1
    WTYSUM[i+1,j+1]=WTYSUM[i+1,j+1]+wtmi1j1

    # Interpolation to P-nodes
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xp[1])/dx)+1
    i=trunc((ym[m]-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xm[m]-xp[j]
    dymi=ym[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update properties
    # i;j Node
    GGGPSUM[i,j]=GGGPSUM[i,j]+1/gggtotalm*wtmij
    SXXSUM[i,j]=SXXSUM[i,j]+sxxm[m]*wtmij
    RHOSUM[i,j]=RHOSUM[i,j]+rhototalm*wtmij
    RHOCPSUM[i,j]=RHOCPSUM[i,j]+rhocptotalm*wtmij
    ALPHASUM[i,j]=ALPHASUM[i,j]+alphasolidm[tm[m]]*wtmij
    ALPHAFSUM[i,j]=ALPHAFSUM[i,j]+alphafluidm[tm[m]]*wtmij
    HRSUM[i,j]=HRSUM[i,j]+hrtotalm*wtmij
    TKSUM[i,j]=TKSUM[i,j]+tkm[m]*rhocptotalm*wtmij
    PHISUM[i,j]=PHISUM[i,j]+phim[m]*wtmij
    WTPSUM[i,j]=WTPSUM[i,j]+wtmij
    # i+1;j Node
    GGGPSUM[i+1,j]=GGGPSUM[i+1,j]+1/gggtotalm*wtmi1j
    SXXSUM[i+1,j]=SXXSUM[i+1,j]+sxxm[m]*wtmi1j
    RHOSUM[i+1,j]=RHOSUM[i+1,j]+rhototalm*wtmi1j
    RHOCPSUM[i+1,j]=RHOCPSUM[i+1,j]+rhocptotalm*wtmi1j
    ALPHASUM[i+1,j]=ALPHASUM[i+1,j]+alphasolidm[tm[m]]*wtmi1j
    ALPHAFSUM[i+1,j]=ALPHAFSUM[i+1,j]+alphafluidm[tm[m]]*wtmi1j
    HRSUM[i+1,j]=HRSUM[i+1,j]+hrtotalm*wtmi1j
    TKSUM[i+1,j]=TKSUM[i+1,j]+tkm[m]*rhocptotalm*wtmi1j
    PHISUM[i+1,j]=PHISUM[i+1,j]+phim[m]*wtmi1j
    WTPSUM[i+1,j]=WTPSUM[i+1,j]+wtmi1j
    # i;j+1 Node
    GGGPSUM[i,j+1]=GGGPSUM[i,j+1]+1/gggtotalm*wtmij1
    SXXSUM[i,j+1]=SXXSUM[i,j+1]+sxxm[m]*wtmij1
    RHOSUM[i,j+1]=RHOSUM[i,j+1]+rhototalm*wtmij1
    RHOCPSUM[i,j+1]=RHOCPSUM[i,j+1]+rhocptotalm*wtmij1
    ALPHASUM[i,j+1]=ALPHASUM[i,j+1]+alphasolidm[tm[m]]*wtmij1
    ALPHAFSUM[i,j+1]=ALPHAFSUM[i,j+1]+alphafluidm[tm[m]]*wtmij1
    HRSUM[i,j+1]=HRSUM[i,j+1]+hrtotalm*wtmij1
    TKSUM[i,j+1]=TKSUM[i,j+1]+tkm[m]*rhocptotalm*wtmij1
    PHISUM[i,j+1]=PHISUM[i,j+1]+phim[m]*wtmij1
    WTPSUM[i,j+1]=WTPSUM[i,j+1]+wtmij1
    # i+1;j+1 Node
    GGGPSUM[i+1,j+1]=GGGPSUM[i+1,j+1]+1/gggtotalm*wtmi1j1
    SXXSUM[i+1,j+1]=SXXSUM[i+1,j+1]+sxxm[m]*wtmi1j1
    RHOSUM[i+1,j+1]=RHOSUM[i+1,j+1]+rhototalm*wtmi1j1
    RHOCPSUM[i+1,j+1]=RHOCPSUM[i+1,j+1]+rhocptotalm*wtmi1j1
    ALPHASUM[i+1,j+1]=ALPHASUM[i+1,j+1]+alphasolidm[tm[m]]*wtmi1j1
    ALPHAFSUM[i+1,j+1]=ALPHAFSUM[i+1,j+1]+alphafluidm[tm[m]]*wtmi1j1
    HRSUM[i+1,j+1]=HRSUM[i+1,j+1]+hrtotalm*wtmi1j1
    TKSUM[i+1,j+1]=TKSUM[i+1,j+1]+tkm[m]*rhocptotalm*wtmi1j1
    PHISUM[i+1,j+1]=PHISUM[i+1,j+1]+phim[m]*wtmi1j1
    WTPSUM[i+1,j+1]=WTPSUM[i+1,j+1]+wtmi1j1
end
# Compute physical properties
# Basic nodes
YNY=zeros(Ny,Nx)
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM[i,j]>0)
            ETA0[i,j]=ETA0SUM[i,j]/WTSUM[i,j]
            ETA[i,j]=ETASUM[i,j]/WTSUM[i,j]
            if(ETA[i,j]<ETA0[i,j])
                YNY[i,j]=1
            end
            GGG[i,j]=1/(GGGSUM[i,j]/WTSUM[i,j])
            SXY0[i,j]=SXYSUM[i,j]/WTSUM[i,j]
            COH[i,j]=COHSUM[i,j]/WTSUM[i,j]
            TEN[i,j]=TENSUM[i,j]/WTSUM[i,j]
            FRI[i,j]=FRISUM[i,j]/WTSUM[i,j]
        end
    end
end
# Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM[i,j]>0)
            RHOX[i,j]=RHOXSUM[i,j]/WTXSUM[i,j]
            RHOFX[i,j]=RHOFXSUM[i,j]/WTXSUM[i,j]
            KX[i,j]=KXSUM[i,j]/WTXSUM[i,j]
            PHIX[i,j]=PHIXSUM[i,j]/WTXSUM[i,j]
            RX[i,j]=RXSUM[i,j]/WTXSUM[i,j]
        end
    end
end
# Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM[i,j]>0)
            RHOY[i,j]=RHOYSUM[i,j]/WTYSUM[i,j]
            RHOFY[i,j]=RHOFYSUM[i,j]/WTYSUM[i,j]
            KY[i,j]=KYSUM[i,j]/WTYSUM[i,j]
            PHIY[i,j]=PHIYSUM[i,j]/WTYSUM[i,j]
            RY[i,j]=RYSUM[i,j]/WTYSUM[i,j]
        end
    end
end
# P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM[i,j]>0)
            GGGP[i,j]=1/(GGGPSUM[i,j]/WTPSUM[i,j])
            SXX0[i,j]=SXXSUM[i,j]/WTPSUM[i,j]
            RHO[i,j]=RHOSUM[i,j]/WTPSUM[i,j]
            RHOCP[i,j]=RHOCPSUM[i,j]/WTPSUM[i,j]
            ALPHA[i,j]=ALPHASUM[i,j]/WTPSUM[i,j]
            ALPHAF[i,j]=ALPHAFSUM[i,j]/WTPSUM[i,j]
            HR[i,j]=HRSUM[i,j]/WTPSUM[i,j]
            PHI[i,j]=PHISUM[i,j]/WTPSUM[i,j]
            BETTAPHI[i,j]=1/GGGP[i,j]*PHI[i,j]
            tk1[i,j]=TKSUM[i,j]/RHOCPSUM[i,j]
        end
    end
end
# Applying thermal boundary conditions for interpolated temperature
# Upper boundary 
tk1[1,2:Nx]=tk1[2,2:Nx]; # Insulating boundary
# Lower boundary 
tk1[Ny1,2:Nx]=tk1[Ny,2:Nx]; # Insulating boundary
# Left boundary
tk1[:,1]=tk1[:,2]; # Insulating boundary
# Right boundary
tk1[:,Nx1]=tk1[:,Nx]; # Insulating boundary



# Gravity solution
# Composing global matrixes LT[], RT[]
# Going through all points of the 2D grid &
# composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        # Define global index in algebraic space
        gk=(j-1)*Ny1+i
        # Distance from the model centre
        rnode=((xp[j]-xsize/2)^2+(yp[i]-ysize/2)^2)^0.5
        # External points
        if(rnode>xsize/2 || i==1 || i==Ny1 || j==1 || j==Nx1)
            # Boundary Condition
            # PHI=0
            LP[gk,gk]=1; # Left part
            RP[gk]=0; # Right part
        else()
        # Internal points: Temperature eq.
        # d2PHI/dx^2+d2PHI/dy^2=2/3*4*G*pi*RHO
        #          PHI2
        #           |
        #           |
        #  PHI1----PHI3----PHI5
        #           |
        #           |
        #          PHI4
        #
        # Density gradients
        dRHOdx=(RHO[i,j+1]-RHO[i,j-1])/2/dx
        dRHOdy=(RHO[i+1,j]-RHO[i-1,j])/2/dy
        # Left part
        LP[gk,gk-Ny1]=1/dx^2; # PHI1
        LP[gk,gk-1]=1/dy^2; # PHI2
        LP[gk,gk]=-2/dx^2-2/dy^2; # PHI3
        LP[gk,gk+1]=1/dy^2; # PHI4
        LP[gk,gk+Ny1]=1/dx^2; # PHI5
        # Right part
        RP[gk]=2/3*4*G*pi*RHO[i,j]
        end
    end
end

# Solving matrixes
SP=LP\RP; # Obtaining algebraic vector of solutions SP[]

# Reload solutions SP[] to geometrical array PHI[]
# Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        # Compute global index
        gk=(j-1)*Ny1+i
        # Reload solution
        FI[i,j]=SP[gk]
    end
end
# Compute gravity acceleration
# gx
for j=1:1:Nx
    for i=1:1:Ny1
        # gx=-dPHI/dx
        gx[i,j]=-(FI[i,j+1]-FI[i,j])/dx
    end
end
# gy
for j=1:1:Nx1
    for i=1:1:Ny
        # gy=-dPHI/dy
        gy[i,j]=-(FI[i+1,j]-FI[i,j])/dy
    end
end





# Try to increase computational Timestep
dt=min(dt*dtkoefup,dtelastic)

# # Set initial viscoplastic viscosity
# if(timestep==1)
#     ETA=ETA0
# end

# Save initial viscoplastic viscosity
ETA00=ETA
# Save initial yielding nodes
YNY00=YNY

# Start Plastic iterations on Nodes until 'End Plastic iterations on Nodes'
if (timestep==1)
    BETTAPHI=zeros(Ny1,Nx1); # No elastic compaction for the first timestep
end
for iplast=1:1:nplast

# Recompute viscosity at pressure nodes
for i=2:1:Ny
    for j=2:1:Nx
        ETAP[i,j]=1/((1/ETA[i-1,j-1]+1/ETA[i,j-1]+1/ETA[i-1,j]+1/ETA[i,j])/4)
        ETAPHI[i,j]=etaphikoef*ETAP[i,j]/PHI[i,j]
    end
end
# Hydro-Mechanical Solution
# Composing global matrixes L[], R[] for Stokes & continuity equations
for j=1:1:Nx1
    for i=1:1:Ny1
        # Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; # Vx solid
        kvy=kvx+1; # Vy solid
        kpm=kvx+2; # Ptotal
        kqx=kvx+3; # qx Darcy
        kqy=kvx+4; # qy Darcy
        kpf=kvx+5; # P fluid
        
        # Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            # Boundary Condition 
            # Ghost unknowns 1*Vx=0
            if(j==Nx1)
                L[kvx,kvx]=1; # Left part
                R[kvx]=0; # Right part
            end
            # Left Boundary
            if(j==1)
                L[kvx,kvx]=1; # Left part
                R[kvx]=vxleft; # Right part
            end
            # Right Boundary
            if(j==Nx)
                L[kvx,kvx]=1; # Left part
                R[kvx]=vxright; # Right part
            end
            # Top boundary
            if(i==1 && j>1 && j<Nx)
                L[kvx,kvx]=1; # Left part
                L[kvx,kvx+6]=bctop; # Left part
                R[kvx]=0; # Right part
            end
            # Top boundary
            if(i==Ny1 && j>1 && j<Nx)
                L[kvx,kvx]=1; # Left part
                L[kvx,kvx-6]=bcbottom; # Left part
                R[kvx]=0; # Right part
            end
        else()
        # Internal points: x-Stokes eq.
        #            Vx2
        #             |
        #        Vy1  |  Vy3
        #             |
        #     Vx1-P1-Vx3-P2-Vx5
        #             |
        #        Vy2  |  Vy4
        #             |
        #            Vx4
        #
        # Computational viscosity
        ETA1=ETA[i-1,j]*GGG[i-1,j]*dt/(GGG[i-1,j]*dt+ETA[i-1,j])
        ETA2=ETA[i,j]*GGG[i,j]*dt/(GGG[i,j]*dt+ETA[i,j])
        ETAP1=ETAP[i,j]*GGGP[i,j]*dt/(GGGP[i,j]*dt+ETAP[i,j])
        ETAP2=ETAP[i,j+1]*GGGP[i,j+1]*dt/(GGGP[i,j+1]*dt+ETAP[i,j+1])
        # Old stresses
        SXY1=SXY0[i-1,j]*ETA[i-1,j]/(GGG[i-1,j]*dt+ETA[i-1,j])
        SXY2=SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dt+ETA[i,j])
        SXX1=SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dt+ETAP[i,j])
        SXX2=SXX0[i,j+1]*ETAP[i,j+1]/(GGGP[i,j+1]*dt+ETAP[i,j+1])
        # Density gradients
        dRHOdx=(RHOX[i,j+1]-RHOX[i,j-1])/2/dx
        dRHOdy=(RHOX[i+1,j]-RHOX[i-1,j])/2/dy
        # Left part
        L[kvx,kvx-Ny1*6]=ETAP1/dx^2; # Vx1
        L[kvx,kvx-6]=ETA1/dy^2; # Vx2
        L[kvx,kvx]=-(ETAP1+ETAP2)/dx^2-  (ETA1+ETA2)/dy^2-  dRHOdx*gx[i,j]*dt; # Vx3
        L[kvx,kvx+6]=ETA2/dy^2; # Vx4
        L[kvx,kvx+Ny1*6]=ETAP2/dx^2; # Vx5
        L[kvx,kvy]=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy2
        L[kvx,kvy+Ny1*6]=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy4
        L[kvx,kvy-6]=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy1
        L[kvx,kvy+Ny1*6-6]=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy3
        L[kvx,kpm]=pscale/dx; # P1
        L[kvx,kpm+Ny1*6]=-pscale/dx; # P2
        # Right part
        R[kvx]=-RHOX[i,j]*gx[i,j]-(SXY2-SXY1)/dy-(SXX2-SXX1)/dx
        end
        
        # Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            # Boundary Condition
            # Ghost unknowns 1*Vx=0
            if(i==Ny1)
                L[kvy,kvy]=1; # Left part
                R[kvy]=0; # Right part
            end
            # Top boundary
            if(i==1)
                L[kvy,kvy]=1; # Left part
                R[kvy]=vytop; # Right part
            end
            # Bottom boundary
            if(i==Ny)
                L[kvy,kvy]=1; # Left part
                R[kvy]=vybottom; # Right part
            end
            # Left boundary
            if(j==1 && i>1 && i<Ny)
                L[kvy,kvy]=1; # Left part
                L[kvy,kvy+6*Ny1]=bcleft; # Left part
                R[kvy]=0; # Right part
            end
            # Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L[kvy,kvy]=1; # Left part
                L[kvy,kvy-6*Ny1]=bcright; # Left part
                R[kvy]=0; # Right part
            end
        else()
        # Internal points: y-Stokes eq.
        #            Vy2
        #             |
        #         Vx1 P1 Vx3
        #             |
        #     Vy1----Vy3----Vy5
        #             |
        #         Vx2 P2 Vx4
        #             |
        #            Vy4
        #
        # Computational viscosity
        ETA1=ETA[i,j-1]*GGG[i,j-1]*dt/(GGG[i,j-1]*dt+ETA[i,j-1])
        ETA2=ETA[i,j]*GGG[i,j]*dt/(GGG[i,j]*dt+ETA[i,j])
        ETAP1=ETAP[i,j]*GGGP[i,j]*dt/(GGGP[i,j]*dt+ETAP[i,j])
        ETAP2=ETAP[i+1,j]*GGGP[i+1,j]*dt/(GGGP[i+1,j]*dt+ETAP[i+1,j])
        # Old stresses
        SXY1=SXY0[i,j-1]*ETA[i,j-1]/(GGG[i,j-1]*dt+ETA[i,j-1])
        SXY2=SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dt+ETA[i,j])
        SYY1=-SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dt+ETAP[i,j])
        SYY2=-SXX0[i+1,j]*ETAP[i+1,j]/(GGGP[i+1,j]*dt+ETAP[i+1,j])
        # Density gradients
        dRHOdx=(RHOY[i,j+1]-RHOY[i,j-1])/2/dx
        dRHOdy=(RHOY[i+1,j]-RHOY[i-1,j])/2/dy
        # Left part
        L[kvy,kvy-Ny1*6]=ETA1/dx^2; # Vy1
        L[kvy,kvy-6]=ETAP1/dy^2; # Vy2
        L[kvy,kvy]=-(ETAP1+ETAP2)/dy^2-  (ETA1+ETA2)/dx^2-  dRHOdy*gy[i,j]*dt; # Vy3
        L[kvy,kvy+6]=ETAP2/dy^2; # Vy4
        L[kvy,kvy+Ny1*6]=ETA2/dx^2; # Vy5
        L[kvy,kvx]=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx3
         L[kvy,kvx+6]=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx4
        L[kvy,kvx-Ny1*6]=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx1
        L[kvy,kvx+6-Ny1*6]=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx2
        L[kvy,kpm]=pscale/dy; # P1
        L[kvy,kpm+6]=-pscale/dy; # P2
        
        # Right part
        R[kvy]=-RHOY[i,j]*gy[i,j]-(SXY2-SXY1)/dx-(SYY2-SYY1)/dy
        end
        
        # P equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1)
            # Boundary Condition
            # 1*P=0
            L[kpm,kpm]=1; # Left part
            R[kpm]=0; # Right part
        else()
        # Internal points: continuity eq.
        # dVx/dx+dVy/dy=0
        #            Vy1
        #             |
        #        Vx1--P--Vx2
        #             |
        #            Vy2
        #
        # Left part
        L[kpm,kvx-Ny1*6]=-1/dx; # Vx1
        L[kpm,kvx]=1/dx; # Vx2
        L[kpm,kvy-6]=-1/dy; # Vy1
        L[kpm,kvy]=1/dy; # Vy2
        L[kpm,kpm]= pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Ptotal
        L[kpm,kpf]=-pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Pfluid
        # Right part
        R[kpm]=(pr0[i,j]-pf0[i,j])/(1-PHI[i,j])*BETTAPHI[i,j]/dt
        end

        # qxDarcy equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            # Boundary Condition
            # 1*qx=0
            L[kqx,kqx]=1; # Left part
            R[kqx]=0; # Right part
            # Top boundary
            if(i==1 && j>1 && j<Nx)
                L[kqx,kqx+6]=bcftop; # Left part
            end
            # Bottom boundary
            if(i==Ny1 && j>1 && j<Nx)
                L[kqx,kqx-6]=bcfbottom; # Left part
            end
        else()
        # Internal points: x-Darcy eq.
        # Rx*qxDarcy+dP/dx=RHOfluid*gx
        #     P1-qxD-P2
        # Left part
        L[kqx,kqx]=RX[i,j]; # qxD
        L[kqx,kpf]=-pscale/dx; # P1
        L[kqx,kpf+Ny1*6]=pscale/dx; # P2
        # Right part
        R[kqx]=RHOFX[i,j]*gx[i,j]
        end
        
        # qyDarcy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            # Boundary Condition
            # 1*Vy=0
            L[kqy,kqy]=1; # Left part
            R[kqy]=0; # Right part
            # Left boundary
            if(j==1 && i>1 && i<Ny)
                L[kqy,kqy+6*Ny1]=bcfleft; # Left part
            end
            # Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L[kqy,kqy-6*Ny1]=bcfright; # Left part
            end
        else()
        # Internal points: y-Stokes eq.
        # Internal points: x-Darcy eq.
        # Rx*qxDarcy+dP/dx=RHOfluid*gx
        #      P1
        #      |
        #     qxD
        #      |
        #      P2
        # Left part
        L[kqy,kqy]=RY[i,j]; # qxD
        L[kqy,kpf]=-pscale/dy; # P1
        L[kqy,kpf+6]=pscale/dy; # P
        # Right part
        R[kqy]=RHOFY[i,j]*gy[i,j]
        end
        
        # Pfluid equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1 || (i==2 && j==2))
            # Boundary Condition
            # 1*Pfluid=0
            L[kpf,kpf]=1; # Left part
            R[kpf]=0; # Right part
            # Real BC
            if(i==2 && j==2)
                L[kpf,kpf]=1*pscale; #Left part
                R[kpf]=psurface; # Right part
            end
        else()
        # Internal points: continuity eq.
        # dqxD/dx+dqyD/dy-(Ptotal-Pfluid)/ETHAphi=0
        #            qyD1
        #              |
        #        qxD1--P--qxD2
        #              |
        #            qyD2
        #
        # Left part
        L[kpf,kqx-Ny1*6]=-1/dx; # qxD1
        L[kpf,kqx]=1/dx; # qxD2
        L[kpf,kqy-6]=-1/dy; # qyD1
        L[kpf,kqy]=1/dy; # qyD2
        L[kpf,kpm]=-pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Ptotal
        L[kpf,kpf]= pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Pfluid
        # Right part
        R[kpf]=-(pr0[i,j]-pf0[i,j])/(1-PHI[i,j])*BETTAPHI[i,j]/dt
        end

    end
end

# 4) Solving matrixes; reloading solution
S=L\R; # Obtaining algebraic vector of solutions S[]
# Reload solutions S[] to vx[], vy[], p[]
# Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        # Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; # Vx solid
        kvy=kvx+1; # Vy solid
        kpm=kvx+2; # Ptotal
        kqx=kvx+3; # qx Darcy
        kqy=kvx+4; # qy Darcy
        kpf=kvx+5; # P fluid
        # Reload solution
        vx[i,j]=S[kvx]
        vy[i,j]=S[kvy]
        pr[i,j]=S[kpm]*pscale
        qxD[i,j]=S[kqx]
        qyD[i,j]=S[kqy]
        pf[i,j]=S[kpf]*pscale
    end
end

# Compute Dln[(1-PHI)/PHI]/Dt
APHI=zeros(Ny1,Nx1)
aphimax=0
for j=2:1:Nx
    for i=2:1:Ny
        APHI[i,j]=((pr[i,j]-pf[i,j])/ETAPHI[i,j]+  ((pr[i,j]-pr0[i,j])-(pf[i,j]-pf0[i,j]))/dt*BETTAPHI[i,j])/(1-PHI[i,j])/PHI[i,j]
        aphimax=max(aphimax,abs(APHI[i,j]))
    end
end

# Compute fluid velocity
# Vx fluid
for j=1:1:Nx
    for i=2:1:Ny
        vxf[i,j]=qxD[i,j]/PHIX[i,j]
    end
end
# Apply BC
# Top
vxf[1,:]=-bcftop*vxf[2,:];    
# Bottom
vxf[Ny1,:]=-bcfbottom*vxf[Ny,:];    
# Vy fluid
for j=2:1:Nx
    for i=1:1:Ny
        vyf[i,j]=qyD[i,j]/PHIY[i,j]
    end
end
# Apply BC
# Left
vyf[:,1]=-bcfleft*vyf[:,2];    
# Right
vyf[:,Nx1]=-bcfright*vyf[:,Nx];     
# Add solid velocity
vxf0=vxf; vxf=vxf+vx
vyf0=vyf; vyf=vyf+vy


# Define displacement timestep dtm
dtm=dt
maxvx=max(max(abs(vx)))
maxvy=max(max(abs(vy)))
if(dtm*maxvx>dxymax*dx)
    dtm=dxymax*dx/maxvx
end
if(dtm*maxvy>dxymax*dy)
    dtm=dxymax*dy/maxvy
end
# Fluid velocity
maxvxf=max(max(abs(vxf)))
maxvyf=max(max(abs(vyf)))
if(dtm*maxvxf>dxymax*dx)
    dtm=dxymax*dx/maxvxf
end
if(dtm*maxvyf>dxymax*dy)
    dtm=dxymax*dy/maxvyf
end
# Porosity change
if(aphimax*dtm>dphimax)
    dtm=dphimax/aphimax
end



# Compute Stress; stress change & strain rate components
# Compute EPSILONxy; SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); # Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); # Stress SIGMAxy, Pa
DSXY=zeros(Ny,Nx); # Stress change SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        # EXY;SXY; DSXY
        EXY[i,j]=0.5*((vx[i+1,j]-vx[i,j])/dy+  (vy[i,j+1]-vy[i,j])/dx)
        SXY[i,j]=2*ETA[i,j]*EXY[i,j]*GGG[i,j]*dtm/(GGG[i,j]*dtm+ETA[i,j])+  SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dtm+ETA[i,j])
        DSXY[i,j]=SXY[i,j]-SXY0[i,j]
    end
end
# Compute EPSILONxx; SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); # Strain rate EPSILONxx, 1/s
EII=zeros(Ny1,Nx1); # Second strain rate invariant, 1/s
SXX=zeros(Ny1,Nx1); # Stress SIGMA'xx, Pa
SII=zeros(Ny1,Nx1); # Second stress invariant, Pa
DSXX=zeros(Ny1,Nx1); # Stress change SIGMA'xx, Pa
DIVV=zeros(Ny1,Nx1); # div[v]
for j=2:1:Nx
    for i=2:1:Ny
        # DIVV
        DIVV[i,j]=(vx[i,j]-vx[i,j-1])/dx+(vy[i,j]-vy[i-1,j])/dy
        # EXX
        EXX[i,j]=((vx[i,j]-vx[i,j-1])/dx-(vy[i,j]-vy[i-1,j])/dy)/2
        # SXX
        SXX[i,j]=2*ETAP[i,j]*EXX[i,j]*GGGP[i,j]*dtm/(GGGP[i,j]*dtm+ETAP[i,j])+  SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dtm+ETAP[i,j])
        DSXX[i,j]=SXX[i,j]-SXX0[i,j]
        # EII
        EII[i,j]=(EXX[i,j]^2+((EXY[i,j]+EXY[i-1,j]+  EXY[i,j-1]+EXY[i-1,j-1])/4)^2)^0.5
        # SII
        SII[i,j]=(SXX[i,j]^2+((SXY[i,j]+SXY[i-1,j]+ SXY[i,j-1]+SXY[i-1,j-1])/4)^2)^0.5
    end
end

# Recompute Dln[(1-PHI)/PHI]/Dt
APHI=zeros(Ny1,Nx1)
for j=2:1:Nx
    for i=2:1:Ny
        APHI[i,j]=((pr[i,j]-pf[i,j])/ETAPHI[i,j]+ ((pr[i,j]-pr0[i,j])-(pf[i,j]-pf0[i,j]))/dt*BETTAPHI[i,j])/(1-PHI[i,j])/PHI[i,j]
    end
end


# Apply Symmetry to Pressure nodes
# External P-nodes: symmetry
# Top
SXX[1,2:Nx]=SXX[2,2:Nx]
APHI[1,2:Nx]=APHI[2,2:Nx];    
PHI[1,2:Nx]=PHI[2,2:Nx];    
pr[1,2:Nx]=pr[2,2:Nx];    
pf[1,2:Nx]=pf[2,2:Nx];    
# Bottom
SXX[Ny1,2:Nx]=SXX[Ny,2:Nx]
APHI[Ny1,2:Nx]=APHI[Ny,2:Nx];    
PHI[Ny1,2:Nx]=PHI[Ny,2:Nx];    
pr[Ny1,2:Nx]=pr[Ny,2:Nx];    
pf[Ny1,2:Nx]=pf[Ny,2:Nx];    
# Left
SXX[:,1]=SXX[:,2]
APHI[:,1]=APHI[:,2];    
PHI[:,1]=PHI[:,2];    
pr[:,1]=pr[:,2];    
pf[:,1]=pf[:,2];    
# Right
SXX[:,Nx1]=SXX[:,Nx]
APHI[:,Nx1]=APHI[:,Nx];    
PHI[:,Nx1]=PHI[:,Nx];    
pr[:,Nx1]=pr[:,Nx];    
pf[:,Nx1]=pf[:,Nx]; 

# Compute solid pressure
ps=(pr-pf.*PHI)./(1-PHI)


# Save nodal stress changes
DSXX0=DSXX
DSXY0=DSXY


# Nodal adjusment
# Update viscosity for yielding
# Basic nodes
ETA5=ETA0
YNY5=zeros(Ny,Nx)
DSY=zeros(Ny,Nx)
ynpl=0
ddd=0
for i=1:1:Ny
  for j=1:1:Nx
    # Compute second stress invariant
    # SXX; pt are averaged from four surrounding pressure nodes
#     SIIB=(SXY[i,j]^2+(SXX[i,j]^2+SXX[i+1,j]^2+SXX[i,j+1]^2+SXX[i+1,j+1]^2)/4)^0.5
    SIIB=(SXY[i,j]^2+((SXX[i,j]+SXX[i+1,j]+SXX[i,j+1]+SXX[i+1,j+1])/4)^2)^0.5
    # Compute second invariant for a purely elastic stress build-up
    siiel=SIIB*(GGG[i,j]*dt+ETA[i,j])/ETA[i,j]
    # Compute total & fluid pressure
    prB=(pr[i,j]+pr[i+1,j]+pr[i,j+1]+pr[i+1,j+1])/4
    pfB=(pf[i,j]+pf[i+1,j]+pf[i,j+1]+pf[i+1,j+1])/4
    # Compute yielding stress
    syieldc=COH[i,j]+FRI[i,j]*(prB-pfB); # Confined fracture
    syieldt=TEN[i,j]+(prB-pfB); # Tensile fracture
    syield=max(minimum(syieldt,syieldc),0); # Non-negative strength requirement
    # Update error for old yielding nodes
    ynn=0
    if(YNY[i,j]>0)
        ynn=1
        DSY[i,j]=SIIB-syield
        ddd=ddd+DSY[i,j]^2
        ynpl=ynpl+1
    end
    # Correcting viscosity for yielding
    if(syield<siiel)
        # New viscosity for the basic node
        etapl=dt*GGG[i,j]*syield/(siiel-syield)
        if(etapl<ETA0[i,j])
            # Recompute nodal visocity
            ETA5[i,j]=etapl^(1-etawt)*ETA[i,j]^etawt
            # Mark yielding nodes
            YNY5[i,j]=1
            # Apply viscosity cutoff values
            if(ETA5[i,j]<etamin)
                ETA5[i,j]=etamin
            elseif[ETA5[i,j]>etamax]
                ETA5[i,j]=etamax
            end
            # Update Error for new yielding nodes
            if(ynn==0)
                DSY[i,j]=SIIB-syield
                ddd=ddd+DSY[i,j]^2
                ynpl=ynpl+1
            end
        else()
            ETA5[i,j]=ETA0[i,j]
            YNY5[i,j]=0
        end
    else()
        ETA5[i,j]=ETA0[i,j]
        YNY5[i,j]=0
    end
  end
end
# Compute yielding error for markers
if(ynpl>0)
    YERRNOD[iplast]=(ddd/ynpl)^0.5
end


# Stop iteration
if(ynpl==0 || iplast==nplast || YERRNOD[iplast]<yerrmax)
    break
# Repeat iteration
else()
    # Decrease computational timestep if too many iterations
    if(trunc(iplast/dtstep)*dtstep==iplast)
        # Decrease timestep
        dt=dt/dtkoef
        # Reset old viscoplastic viscosity
        ETA=ETA00
        YNY=YNY00
    else()
        # Use new viscoplastic viscosity
        ETA=ETA5
        YNY=YNY5
    end
end
# End Plastic iterations on Nodes
end



# Interpolate new viscoplastic viscosity to markers
for m=1:1:marknum
    # Interpolation viscosity from basic nodes
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-x[1])/dx)+1
    i=trunc((ym[m]-y[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<1)
        i=1
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-x[j]
    dymi=ym[m]-y[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Matrix viscosity
    if(tm[m]<3)
        # Rocks
        etasolidcur=etasolidm[tm[m]]
        if(tkm[m]>tmsilicate)
            etasolidcur=etasolidmm[tm[m]]
        end
        etatotalm=etasolidcur;#*exp(-28*phim[m])
    else()
        # Sticky air
        etatotalm=etasolidm[tm[m]]
    end
    if(YNY[i,j]>0 || YNY[i+1,j]>0 || YNY[i,j+1]>0 || YNY[i+1,j+1]>0)
#         etavpm[m]=ETA[i,j]*wtmij+ETA[i+1,j]*wtmi1j+...
#                 ETA[i,j+1]*wtmij1+ETA[i+1,j+1]*wtmi1j1
#         etavpm[m]=1/(1/ETA[i,j]*wtmij+1/ETA[i+1,j]*wtmi1j+...
#                 1/ETA[i,j+1]*wtmij1+1/ETA[i+1,j+1]*wtmi1j1)
        etavpm[m]=1/(YNY[i,j]/ETA[i,j]*wtmij+YNY[i+1,j]/ETA[i+1,j]*wtmi1j+ YNY[i,j+1]/ETA[i,j+1]*wtmij1+YNY[i+1,j+1]/ETA[i+1,j+1]*wtmi1j1)
        if(etavpm[m]>=etatotalm)
            etavpm[m]=etatotalm
        end
    else()
        etavpm[m]=etatotalm
    end
end



# Apply subgrid stress diffusion to markers
if(dsubgrids>0)
SXYSUM=zeros(Ny,Nx)
WTSUM=zeros(Ny,Nx)
SXXSUM=zeros(Ny1,Nx1)
WTPSUM=zeros(Ny1,Nx1)
for m=1:1:marknum
    # SIGMA'xx
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xp[1])/dx)+1
    i=trunc((ym[m]-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xm[m]-xp[j]
    dymi=ym[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute marker-node SIGMA'xx difference
    dsxxm0=sxxm[m]-(SXX0[i,j]*wtmij+SXX0[i+1,j]*wtmi1j+ SXX0[i,j+1]*wtmij1+SXX0[i+1,j+1]*wtmi1j1)
    # Relax stress difference
    dsxxm1=dsxxm0*exp(-dsubgrids*dtm/(etam[tm[m]]/gggm[tm[m]]))
    # Correct marker stress
    ddsxxm=dsxxm1-dsxxm0
    sxxm[m]=sxxm[m]+ddsxxm
    # Update subgrid diffusion on nodes
    # i;j Node
    SXXSUM[i,j]=SXXSUM[i,j]+ddsxxm*wtmij
    WTPSUM[i,j]=WTPSUM[i,j]+wtmij
    # i+1;j Node
    SXXSUM[i+1,j]=SXXSUM[i+1,j]+ddsxxm*wtmi1j
    WTPSUM[i+1,j]=WTPSUM[i+1,j]+wtmi1j
    # i;j+1 Node
    SXXSUM[i,j+1]=SXXSUM[i,j+1]+ddsxxm*wtmij1
    WTPSUM[i,j+1]=WTPSUM[i,j+1]+wtmij1
    # i+1;j+1 Node
    SXXSUM[i+1,j+1]=SXXSUM[i+1,j+1]+ddsxxm*wtmi1j1
    WTPSUM[i+1,j+1]=WTPSUM[i+1,j+1]+wtmi1j1

    # SIGMAxy
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-x[1])/dx)+1
    i=trunc((ym[m]-y[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<1)
        i=1
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-x[j]
    dymi=ym[m]-y[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute marker-node SIGMAxy difference
    dsxym0=sxym[m]-(SXY0[i,j]*wtmij+SXY0[i+1,j]*wtmi1j+  SXY0[i,j+1]*wtmij1+SXY0[i+1,j+1]*wtmi1j1)
    # Relax stress difference
    dsxym1=dsxym0*exp(-dsubgrids*dtm/(etam[tm[m]]/gggm[tm[m]]))
    # Correct marker stress
    ddsxym=dsxym1-dsxym0
    sxym[m]=sxym[m]+ddsxym
    # Update subgrid diffusion on nodes
    # i;j Node
    SXYSUM[i,j]=SXYSUM[i,j]+ddsxym*wtmij
    WTSUM[i,j]=WTSUM[i,j]+wtmij
    # i+1;j Node
    SXYSUM[i+1,j]=SXYSUM[i+1,j]+ddsxym*wtmi1j
    WTSUM[i+1,j]=WTSUM[i+1,j]+wtmi1j
    # i;j+1 Node
    SXYSUM[i,j+1]=SXYSUM[i,j+1]+ddsxym*wtmij1
    WTSUM[i,j+1]=WTSUM[i,j+1]+wtmij1
    # i+1;j+1 Node
    SXYSUM[i+1,j+1]=SXYSUM[i+1,j+1]+ddsxym*wtmi1j1
    WTSUM[i+1,j+1]=WTSUM[i+1,j+1]+wtmi1j1
end
# Compute DSXXsubgrid
DSXXsubgrid=zeros(Ny1,Nx1)
# P-nodes
for j=2:1:Nx
    for i=2:1:Ny
        if(WTPSUM[i,j]>0)
            DSXXsubgrid[i,j]=SXXSUM[i,j]/WTPSUM[i,j]
        end
    end
end
# Correct DSXX
DSXX=DSXX-DSXXsubgrid
# Compute DSXYsubgrid
DSXYsubgrid=zeros(Ny,Nx)
# Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM[i,j]>0)
            DSXYsubgrid[i,j]=SXYSUM[i,j]/WTSUM[i,j]
        end
    end
end
# Correct DSXY
DSXY=DSXY-DSXYsubgrid
end

# Interpolate DSXX; DSXY to markers
for m=1:1:marknum
    # SIGMA'xx
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xp[1])/dx)+1
    i=trunc((ym[m]-yp[1])/dy)+1
    if(j<2)
        j=2
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<2)
        i=2
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-xp[j]
    dymi=ym[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update marker by SIGMA'xx change 
    sxxm[m]=sxxm[m]+(DSXX[i,j]*wtmij+DSXX[i+1,j]*wtmi1j+ DSXX[i,j+1]*wtmij1+DSXX[i+1,j+1]*wtmi1j1)

    # SIGMAxy
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-x[1])/dx)+1
    i=trunc((ym[m]-y[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<1)
        i=1
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-x[j]
    dymi=ym[m]-y[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update marker by SIGMA'xx change 
    sxym[m]=sxym[m]+(DSXY[i,j]*wtmij+DSXY[i+1,j]*wtmi1j+ DSXY[i,j+1]*wtmij1+DSXY[i+1,j+1]*wtmi1j1)
end


# Compute shear heating HS in Temperature/Pressure nodes
HS=zeros(Ny1,Nx1); # Adiabatic heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        # Average SXY*EXY
        SXYEXY=(SXY[i,j]^2/ETA[i,j]+SXY[i-1,j]^2/ETA[i-1,j]+ SXY[i,j-1]^2/ETA[i,j-1]+SXY[i-1,j-1]^2/ETA[i-1,j-1])/4
        # HS
        HS[i,j]=SXX[i,j]^2/ETAP[i,j]+SXYEXY+ (pr[i,j]-pf[i,j])^2/(1-PHI[i,j])/ETAPHI[i,j]+ (RX[i,j-1]*qxD[i,j-1]^2+RX[i,j]*qxD[i,j]^2)/2+ (RY[i-1,j]*qyD[i-1,j]^2+RY[i,j]*qyD[i,j]^2)/2
    end
end

# Compute adiabatic heating HA in Temperature/Pressure nodes
# Compute solid pressure
if (timestep==1)
    pr0=pr; # No total pressure change for the first timestep
    pf0=pf; # No fluid pressure change for the first timestep
    ps0=pf; # No solid pressure change for the first timestep
end
# Old solid pressure
HA=zeros(Ny1,Nx1); # Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        # HA
        # Indirect calculation of dpdt
        # Average vy; vx; vxf; vyf
        VXP=(vx[i,j]+vx[i,j-1])/2
        VYP=(vy[i,j]+vy[i-1,j])/2
        VXFP=(vxf[i,j]+vxf[i,j-1])/2
        VYFP=(vyf[i,j]+vyf[i-1,j])/2
        # Evaluate DPsolid/Dt with upwind differences
        if(VXP<0)
            dpsdx=(ps[i,j]-ps[i,j-1])/dx
        else()
            dpsdx=(ps[i,j+1]-ps[i,j])/dx
        end
        if(VYP<0)
            dpsdy=(ps[i,j]-ps[i-1,j])/dy
        else()
            dpsdy=(ps[i+1,j]-ps[i,j])/dy
        end
        dpsdt=VXP*dpsdx+VYP*dpsdy
        # Evaluate DPfluid/Dt with upwind differences
        if(VXFP>0)
            dpfdx=(pf[i,j]-pf[i,j-1])/dx
        else()
            dpfdx=(pf[i,j+1]-pf[i,j])/dx
        end
        if(VYFP>0)
            dpfdy=(pf[i,j]-pf[i-1,j])/dy
        else()
            dpfdy=(pf[i+1,j]-pf[i,j])/dy
        end
        dpfdt=VXFP*dpsdx+VYFP*dpsdy
#         # Direct calculation of dpdt
#         dpsdt=(ps[i,j]-ps0[i,j])/dt
#         dpfdt=(pf[i,j]-pf0[i,j])/dt
        # HA
        HA[i,j]=(1-PHI[i,j])*tk1[i,j]*ALPHA[i,j]*dpsdt+ PHI[i,j]*tk1[i,j]*ALPHAF[i,j]*dpfdt
    end
end


# Thermal iterations
tk0=tk1
dtt=dtm
dttsum=0
titer=1
while(dttsum<dtm)
# Composing global matrixes LT[], RT[]
# Going through all points of the 2D grid &
# composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        # Define global index in algebraic space
        gk=(j-1)*Ny1+i
        # External points
        if(i==1 || i==Ny1 || j==1 || j==Nx1)
            # Boundary Condition
            # Top BC: T=273
            if(i==1 && j>1 && j<Nx1)
                LT[gk,gk]=1; # Left part
                LT[gk,gk+1]=-1; # Left part
                RT[gk]=0; # Right part
            end
            # Bottom BC: T=1500
            if(i==Ny1 && j>1 && j<Nx1)
                LT[gk,gk]=1; # Left part
                LT[gk,gk-1]=-1; # Left part
                RT[gk]=0; # Right part
            end
            # Left BC: dT/dx=0
            if(j==1)
                LT[gk,gk]=1; # Left part
                LT[gk,gk+Ny1]=-1; # Left part
                RT[gk]=0; # Right part
            end
            # Right BC: dT/dx=0
            if(j==Nx1)
                LT[gk,gk]=1; # Left part
                LT[gk,gk-Ny1]=-1; # Left part
                RT[gk]=0; # Right part
            end
        else()
        # Internal points: Temperature eq.
        # RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
        #          Tdt2
        #           |
        #          Ky1
        #           |
        #Tdt1-Kx1-T03;Tdt3-Kx2-Tdt5
        #           |
        #          Ky2
        #           |
        #          Tdt4
        #
        # Left part
        Kx1=KX[i,j-1]; 
        Kx2=KX[i,j]; 
        Ky1=KY[i-1,j]; 
        Ky2=KY[i,j]; 
        LT[gk,gk-Ny1]=-Kx1/dx^2; # T1
        LT[gk,gk-1]=-Ky1/dy^2; # FI2
        LT[gk,gk]=RHOCP[i,j]/dtt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; # FI3
        LT[gk,gk+1]=-Ky2/dy^2; # FI4
        LT[gk,gk+Ny1]=-Kx2/dx^2; # FI5
        # Right part
        RT[gk]=RHOCP[i,j]/dtt*tk1[i,j]+HR[i,j]+HA[i,j]+HS[i,j]
        end
    end
end

# Solving matrixes
ST=LT\RT; # Obtaining algebraic vector of solutions ST[]

# Reload solutions ST[] to geometrical array Tdt[]
# Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        # Compute global index
        gk=(j-1)*Ny1+i
        # Reload solution
        tk2[i,j]=ST[gk]
    end
end
# Compute DT
DT=tk2-tk1
titer
dtt
if(titer==1)
    # Apply thermal timestepping condition
    maxDTcurrent=max(max(abs(DT)))
    if(maxDTcurrent>DTmax)
        dtt=dtt/maxDTcurrent*DTmax
    else()
        dttsum=dttsum+dtt; # Update dttsum
    end
else()
    dttsum=dttsum+dtt; # Update dttsum
    # Adjust timestep
    if(dtt>dtm-dttsum)
        dtt=dtm-dttsum
    end
end

dttsum

titer=titer+1; # Update iteration counter
end
# Compute/save overall temperature changes
DT=tk2-tk0
DT0=DT


# Apply subgrid temperature diffusion on markers
if(dsubgridt>0)
TKSUM=zeros(Ny1,Nx1)
RHOCPSUM=zeros(Ny1,Nx1)
for m=1:1:marknum
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xp[1])/dx)+1
    i=trunc((ym[m]-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xm[m]-xp[j]
    dymi=ym[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute marker-node T difference
    dtkm0=tkm[m]-(tk1[i,j]*wtmij+tk1[i+1,j]*wtmi1j+ tk1[i,j+1]*wtmij1+tk1[i+1,j+1]*wtmi1j1)
    # Compute marker parameters
    if(tm[m]<3)
        # Rocks
        rhocptotalm=rhocpsolidm[tm[m]]*(1-phim[m])+rhocpfluid*phim[m]
        ktotalm=(ksolidm[tm[m]]*kfluid/2+((ksolidm[tm[m]]*(3*phim[m]-2)+ kfluid*(1-3*phim[m]))^2)/16)^0.5-(ksolidm[tm[m]]*(3*phim[m]-2)+ kfluid*(1-3*phim[m]))/4
    else()
        # Sticky air
        rhocptotalm=rhocpsolidm[tm[m]]
        ktotalm=ksolidm[tm[m]]
    end    # Relax temperature difference
    dtkm1=dtkm0*exp(-dsubgridt*ktotalm*dtm/rhocptotalm*(2/dx^2+2/dy^2))
    # Correct marker temperature
    ddtkm=dtkm1-dtkm0
    tkm[m]=tkm[m]+ddtkm
    # Update subgrid diffusion on nodes
    # i;j Node
    TKSUM[i,j]=TKSUM[i,j]+ddtkm*rhocpm[tm[m]]*wtmij
    RHOCPSUM[i,j]=RHOCPSUM[i,j]+rhocpm[tm[m]]*wtmij
    # i+1;j Node
    TKSUM[i+1,j]=TKSUM[i+1,j]+ddtkm*rhocpm[tm[m]]*wtmi1j
    RHOCPSUM[i+1,j]=RHOCPSUM[i+1,j]+rhocpm[tm[m]]*wtmi1j
    # i;j+1 Node
    TKSUM[i,j+1]=TKSUM[i,j+1]+ddtkm*rhocpm[tm[m]]*wtmij1
    RHOCPSUM[i,j+1]=RHOCPSUM[i,j+1]+rhocpm[tm[m]]*wtmij1
    # i+1;j+1 Node
    TKSUM[i+1,j+1]=TKSUM[i+1,j+1]+ddtkm*rhocpm[tm[m]]*wtmi1j1
    RHOCPSUM[i+1,j+1]=RHOCPSUM[i+1,j+1]+rhocpm[tm[m]]*wtmi1j1
end
# Compute DTsubgrid
DTsubgrid=zeros(Ny1,Nx1)
# P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOCPSUM[i,j]>0)
            DTsubgrid[i,j]=TKSUM[i,j]/RHOCPSUM[i,j]
        end
    end
end
# Correct DT
DT=DT-DTsubgrid
end


# Interpolate DT to markers
for m=1:1:marknum
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xp[1])/dx)+1
    i=trunc((ym[m]-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xm[m]-xp[j]
    dymi=ym[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Update properties
    tkm[m]=tkm[m]+DT[i,j]*wtmij+DT[i+1,j]*wtmi1j+ DT[i,j+1]*wtmij1+DT[i+1,j+1]*wtmi1j1
    # Interpolate tk2 at 1st timestep
    if(timestep==1)
        tkm[m]=tk2[i,j]*wtmij+tk2[i+1,j]*wtmi1j+ tk2[i,j+1]*wtmij1+tk2[i+1,j+1]*wtmi1j1
    end
end


    
# Update porosity on markers [not for sticky air]
for m=1:1:marknum
    if(tm[m]<3)
        # Interpolate APHI
        # Define i;j indexes for the upper left node
        j=trunc((xm[m]-xp[1])/dx)+1
        i=trunc((ym[m]-yp[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx]
            j=Nx
        end
        if(i<1)
            i=1
        elseif[i>Ny]
            i=Ny
        end
        # Compute distances
        dxmj=xm[m]-xp[j]
        dymi=ym[m]-yp[i]
        # Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy)
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy)
        wtmi1j1=(dxmj/dx)*(dymi/dy)
        # Compute Dln[(1-PHI)/PHI]/Dt
        aphim=APHI[i,j]*wtmij+APHI[i+1,j]*wtmi1j+ APHI[i,j+1]*wtmij1+APHI[i+1,j+1]*wtmi1j1
        # Change Porosity
        phim[m]=phim[m]/((1-phim[m])*exp(aphim*dtm)+phim[m])
        if(phim[m]<phimin)
            phim[m]=phimin
        elseif[phim[m]>phimax]
            phim[m]=phimax
        end
    end
end

# Compute fluid velocity in pressure nodes
# vxpf
for j=2:1:Nx
    for i=2:1:Ny
        vxpf[i,j]=(vxf[i,j]+vxf[i,j-1])/2
    end
end
# Apply BC
# Top
vxpf[1,2:Nx-1]=-bcftop*vxpf[2,2:Nx-1];    
# Bottom
vxpf[Ny1,2:Nx-1]=-bcfbottom*vxpf[Ny,2:Nx-1];    
# Left
vxpf[:,1]=2*vxleft-vxpf[:,2]
# Right
vxpf[:,Nx1]=2*vxright-vxpf[:,Nx]
# vypf
for j=2:1:Nx
    for i=2:1:Ny
        vypf[i,j]=(vyf[i,j]+vyf[i-1,j])/2
    end
end    
# Apply BC
# Left
vypf[2:Ny-1,1]=-bcfleft*vypf[2:Ny-1,2];    
# Right
vypf[2:Ny-1,Nx1]=-bcfright*vypf[2:Ny-1,Nx]; # Free slip    
# Top
vypf[1,:]=2*vytop-vypf[2,:]
# Bottom
vypf[Ny1,:]=2*vybottom-vypf[Ny,:]

# Compute velocity in pressure nodes
# vx
for j=2:1:Nx
    for i=2:1:Ny
        vxp[i,j]=(vx[i,j]+vx[i,j-1])/2
    end
end
# Apply BC
# Top
vxp[1,2:Nx-1]=-bctop*vxp[2,2:Nx-1];    
# Bottom
vxp[Ny1,2:Nx-1]=-bcbottom*vxp[Ny,2:Nx-1];    
# Left
vxp[:,1]=2*vxleft-vxp[:,2]
# Right
vxp[:,Nx1]=2*vxright-vxp[:,Nx]
# vy
for j=2:1:Nx
    for i=2:1:Ny
        vyp[i,j]=(vy[i,j]+vy[i-1,j])/2
    end
end    
# Apply BC
# Left
vyp[2:Ny-1,1]=-bcleft*vyp[2:Ny-1,2];    
# Right
vyp[2:Ny-1,Nx1]=-bcright*vyp[2:Ny-1,Nx]; # Free slip    
# Top
vyp[1,:]=2*vytop-vyp[2,:]
# Bottom
vyp[Ny1,:]=2*vybottom-vyp[Ny,:]

# Compute rotation rate wyx=1/2[dVy/dx-dVx/dy] for basic nodes
for i=1:1:Ny
    for j=1:1:Nx
        wyx[i,j]=0.5*((vy[i,j+1]-vy[i,j])/dx-(vx[i+1,j]-vx[i,j])/dy)
    end
end

# Move markers with 4th order Runge-Kutta
vxm=zeros(4,1)
vym=zeros(4,1)
for m=1:1:marknum
    
    # Interpolate solid temperature for the initial marker location
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xp[1])/dx)+1
    i=trunc((ym[m]-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xm[m]-xp[j]
    dymi=ym[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute Tsolid
    tksm0=tk2[i,j]*wtmij+tk2[i+1,j]*wtmi1j+ tk2[i,j+1]*wtmij1+tk2[i+1,j+1]*wtmi1j1;    
        
    # Interpolate local rotation rate
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-x[1])/dx)+1
    i=trunc((ym[m]-y[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx-1]
        j=Nx-1
    end
    if(i<1)
        i=1
    elseif[i>Ny-1]
        i=Ny-1
    end
    # Compute distances
    dxmj=xm[m]-x[j]
    dymi=ym[m]-y[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute vx velocity
    omegam=wyx[i,j]*wtmij+wyx[i+1,j]*wtmi1j+ wyx[i,j+1]*wtmij1+wyx[i+1,j+1]*wtmi1j1
    # Analytical stress rotation using SIGMA"xx=-SIGMA"yy
    THETA=dtm*omegam; # Incremental rotation angle()
    sxxmnew=sxxm[m]*cos(THETA)^2-sxxm[m]*sin(THETA)^2-sxym[m]*sin(2*THETA)
    sxymnew=sxxm[m]*sin(2*THETA)+sxym[m]*cos(2*THETA)
    sxxm[m]=sxxmnew; sxym[m]=sxymnew;    
    
    # Save initial marker coordinates
    xA=xm[m]
    yA=ym[m]
    for rk=1:1:4
        # Interpolate vx
        # Define i;j indexes for the upper left node
        j=trunc((xm[m]-xvx[1])/dx)+1
        i=trunc((ym[m]-yvx[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx-1]
            j=Nx-1
        end
        if(i<1)
            i=1
        elseif[i>Ny]
            i=Ny
        end
        # Compute distances
        dxmj=xm[m]-xvx[j]
        dymi=ym[m]-yvx[i]
        # Compute weights
        # Compute vx velocity for the top & bottom of the cell()
        vxm13=vx[i,j]*(1-dxmj/dx)+vx[i,j+1]*dxmj/dx
        vxm24=vx[i+1,j]*(1-dxmj/dx)+vx[i+1,j+1]*dxmj/dx
        # Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx[i,j]-2*vx[i,j+1]+vx[i,j+2])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx[i+1,j]-2*vx[i+1,j+1]+vx[i+1,j+2])
            end
        else()
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx[i,j-1]-2*vx[i,j]+vx[i,j+1])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx[i+1,j-1]-2*vx[i+1,j]+vx[i+1,j+1])
            end
        end
        # Compute vx
        vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24
        
        # Interpolate vy
        # Define i;j indexes for the upper left node
        j=trunc((xm[m]-xvy[1])/dx)+1
        i=trunc((ym[m]-yvy[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx]
            j=Nx
        end
        if(i<1)
            i=1
        elseif[i>Ny-1]
            i=Ny-1
        end
        # Compute distances
        dxmj=xm[m]-xvy[j]
        dymi=ym[m]-yvy[i]
        # Compute weights
        # Compute vy velocity for the left & right of the cell()
        vym12=vy[i,j]*(1-dymi/dy)+vy[i+1,j]*dymi/dy
        vym34=vy[i,j+1]*(1-dymi/dy)+vy[i+1,j+1]*dymi/dy
        # Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy[i,j]-2*vy[i+1,j]+vy[i+2,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy[i,j+1]-2*vy[i+1,j+1]+vy[i+2,j+1])
            end      
        else()
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j]-2*vy[i,j]+vy[i+1,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j+1]-2*vy[i,j+1]+vy[i+1,j+1])
            end
        end
        # Compute vy
        vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34
        
        # Change coordinates to obtain B;C;D points
        if(rk==1 || rk==2)
            xm[m]=xA+dtm/2*vxm[rk]
            ym[m]=yA+dtm/2*vym[rk]
        elseif[rk==3]
            xm[m]=xA+dtm*vxm[rk]
            ym[m]=yA+dtm*vym[rk]
        end
    end
    # Restore initial coordinates
    xm[m]=xA
    ym[m]=yA
    # Compute effective velocity
    vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
    vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
    # Move markers
    xm[m]=xm[m]+dtm*vxmeff
    ym[m]=ym[m]+dtm*vymeff
    
    # Backtracing markers with fluid velocity
    xcur=xm[m]
    ycur=ym[m]
    xA=xcur
    yA=ycur
    for rk=1:1:4
        # Interpolate vx
        # Define i;j indexes for the upper left node
        j=trunc((xcur-xvx[1])/dx)+1
        i=trunc((ycur-yvx[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx-1]
            j=Nx-1
        end
        if(i<1)
            i=1
        elseif[i>Ny]
            i=Ny
        end
        # Compute distances
        dxmj=xcur-xvx[j]
        dymi=ycur-yvx[i]
        # Compute weights
        # Compute vx velocity for the top & bottom of the cell()
        vxm13=vxf[i,j]*(1-dxmj/dx)+vxf[i,j+1]*dxmj/dx
        vxm24=vxf[i+1,j]*(1-dxmj/dx)+vxf[i+1,j+1]*dxmj/dx
        # Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf[i,j]-2*vxf[i,j+1]+vxf[i,j+2])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf[i+1,j]-2*vxf[i+1,j+1]+vxf[i+1,j+2])
            end
        else()
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf[i,j-1]-2*vxf[i,j]+vxf[i,j+1])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf[i+1,j-1]-2*vxf[i+1,j]+vxf[i+1,j+1])
            end
        end
        # Compute vx
        vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24
        
        # Interpolate vy
        # Define i;j indexes for the upper left node
        j=trunc((xcur-xvy[1])/dx)+1
        i=trunc((ycur-yvy[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx]
            j=Nx
        end
        if(i<1)
            i=1
        elseif[i>Ny-1]
            i=Ny-1
        end
        # Compute distances
        dxmj=xcur-xvy[j]
        dymi=ycur-yvy[i]
        # Compute weights
        # Compute vy velocity for the left & right of the cell()
        vym12=vyf[i,j]*(1-dymi/dy)+vyf[i+1,j]*dymi/dy
        vym34=vyf[i,j+1]*(1-dymi/dy)+vyf[i+1,j+1]*dymi/dy
        # Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf[i,j]-2*vyf[i+1,j]+vyf[i+2,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf[i,j+1]-2*vyf[i+1,j+1]+vyf[i+2,j+1])
            end      
        else()
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf[i-1,j]-2*vyf[i,j]+vyf[i+1,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf[i-1,j+1]-2*vyf[i,j+1]+vyf[i+1,j+1])
            end
        end
        # Compute vy
        vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34
        
        # Change coordinates to obtain B;C;D points
        if(rk==1 || rk==2)
            xcur=xA-dtm/2*vxm[rk]
            ycur=yA-dtm/2*vym[rk]
        elseif[rk==3]
            xcur=xA-dtm*vxm[rk]
            ycur=yA-dtm*vym[rk]
        end
    end
    # Compute effective velocity
    vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
    vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
    # Trace the node backward
    xcur=xA-dtm*vxmeff
    ycur=yA-dtm*vymeff
    # Interpolate fluid temperature
    # Define i;j indexes for the upper left node
    j=trunc((xcur-xp[1])/dx)+1
    i=trunc((ycur-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xcur-xp[j]
    dymi=ycur-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute nodal Tfluid
    tkfm0=tk2[i,j]*wtmij+tk2[i+1,j]*wtmi1j+ tk2[i,j+1]*wtmij1+tk2[i+1,j+1]*wtmi1j1
    # Compute Tfluid-Tsolid for the marker
    dtkfsm=tkfm0-tksm0
    # Correct marker temperature
    tkm[m]=((1-phim[m])*tkm[m]*rhocpsolidm[tm[m]]+ phim[m]*(tkm[m]+dtkfsm)*rhocpfluidm[tm[m]])/ ((1-phim[m])*rhocpsolidm[tm[m]]+phim[m]*rhocpfluidm[tm[m]])
end  

# Backtracing Pressure nodes: Ptotal
# Backtracing is based on 4th order Runge-Kutta
vxm=zeros(4,1)
vym=zeros(4,1)
pr0=pr
ps0=ps
for jj=2:1:Nx
for ii=2:1:Ny
    # Save initial nodal coordinates
    xcur=xp[jj]
    ycur=yp[ii]
    xA=xcur
    yA=ycur
    for rk=1:1:4
        # Interpolate vx
        # Define i;j indexes for the upper left node
        j=trunc((xcur-xvx[1])/dx)+1
        i=trunc((ycur-yvx[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx-1]
            j=Nx-1
        end
        if(i<1)
            i=1
        elseif[i>Ny]
            i=Ny
        end
        # Compute distances
        dxmj=xcur-xvx[j]
        dymi=ycur-yvx[i]
        # Compute weights
        # Compute vx velocity for the top & bottom of the cell()
        vxm13=vx[i,j]*(1-dxmj/dx)+vx[i,j+1]*dxmj/dx
        vxm24=vx[i+1,j]*(1-dxmj/dx)+vx[i+1,j+1]*dxmj/dx
        # Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx[i,j]-2*vx[i,j+1]+vx[i,j+2])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx[i+1,j]-2*vx[i+1,j+1]+vx[i+1,j+2])
            end
        else()
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx[i,j-1]-2*vx[i,j]+vx[i,j+1])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx[i+1,j-1]-2*vx[i+1,j]+vx[i+1,j+1])
            end
        end
        # Compute vx
        vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24
        
        # Interpolate vy
        # Define i;j indexes for the upper left node
        j=trunc((xcur-xvy[1])/dx)+1
        i=trunc((ycur-yvy[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx]
            j=Nx
        end
        if(i<1)
            i=1
        elseif[i>Ny-1]
            i=Ny-1
        end
        # Compute distances
        dxmj=xcur-xvy[j]
        dymi=ycur-yvy[i]
        # Compute weights
        # Compute vy velocity for the left & right of the cell()
        vym12=vy[i,j]*(1-dymi/dy)+vy[i+1,j]*dymi/dy
        vym34=vy[i,j+1]*(1-dymi/dy)+vy[i+1,j+1]*dymi/dy
        # Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy[i,j]-2*vy[i+1,j]+vy[i+2,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy[i,j+1]-2*vy[i+1,j+1]+vy[i+2,j+1])
            end      
        else()
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j]-2*vy[i,j]+vy[i+1,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j+1]-2*vy[i,j+1]+vy[i+1,j+1])
            end
        end
        # Compute vy
        vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34
        
        # Change coordinates to obtain B;C;D points
        if(rk==1 || rk==2)
            xcur=xA-dtm/2*vxm[rk]
            ycur=yA-dtm/2*vym[rk]
        elseif[rk==3]
            xcur=xA-dtm*vxm[rk]
            ycur=yA-dtm*vym[rk]
        end
    end
    # Compute effective velocity
    vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
    vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
    # Trace the node backward
    xcur=xA-dtm*vxmeff
    ycur=yA-dtm*vymeff
    # Interpolate nodal property
    # SIGMA'xx; P
    # Define i;j indexes for the upper left node
    j=trunc((xcur-xp[1])/dx)+1
    i=trunc((ycur-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xcur-xp[j]
    dymi=ycur-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute nodal total pressure
    pr0[ii,jj]=pr[i,j]*wtmij+pr[i+1,j]*wtmi1j+ pr[i,j+1]*wtmij1+pr[i+1,j+1]*wtmi1j1
    ps0[ii,jj]=ps[i,j]*wtmij+ps[i+1,j]*wtmi1j+ ps[i,j+1]*wtmij1+ps[i+1,j+1]*wtmi1j1
end
end

# Backtracing Pressure nodes: Pfluid
pf0=pf
for jj=2:1:Nx
for ii=2:1:Ny
    # Save initial nodal coordinates
    xcur=xp[jj]
    ycur=yp[ii]
    xA=xcur
    yA=ycur
    for rk=1:1:4
        # Interpolate vx
        # Define i;j indexes for the upper left node
        j=trunc((xcur-xvx[1])/dx)+1
        i=trunc((ycur-yvx[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx-1]
            j=Nx-1
        end
        if(i<1)
            i=1
        elseif[i>Ny]
            i=Ny
        end
        # Compute distances
        dxmj=xcur-xvx[j]
        dymi=ycur-yvx[i]
        # Compute weights
        # Compute vx velocity for the top & bottom of the cell()
        vxm13=vxf[i,j]*(1-dxmj/dx)+vxf[i,j+1]*dxmj/dx
        vxm24=vxf[i+1,j]*(1-dxmj/dx)+vxf[i+1,j+1]*dxmj/dx
        # Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf[i,j]-2*vxf[i,j+1]+vxf[i,j+2])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf[i+1,j]-2*vxf[i+1,j+1]+vxf[i+1,j+2])
            end
        else()
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf[i,j-1]-2*vxf[i,j]+vxf[i,j+1])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf[i+1,j-1]-2*vxf[i+1,j]+vxf[i+1,j+1])
            end
        end
        # Compute vx
        vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24
        
        # Interpolate vy
        # Define i;j indexes for the upper left node
        j=trunc((xcur-xvy[1])/dx)+1
        i=trunc((ycur-yvy[1])/dy)+1
        if(j<1)
            j=1
        elseif[j>Nx]
            j=Nx
        end
        if(i<1)
            i=1
        elseif[i>Ny-1]
            i=Ny-1
        end
        # Compute distances
        dxmj=xcur-xvy[j]
        dymi=ycur-yvy[i]
        # Compute weights
        # Compute vy velocity for the left & right of the cell()
        vym12=vyf[i,j]*(1-dymi/dy)+vyf[i+1,j]*dymi/dy
        vym34=vyf[i,j+1]*(1-dymi/dy)+vyf[i+1,j+1]*dymi/dy
        # Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf[i,j]-2*vyf[i+1,j]+vyf[i+2,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf[i,j+1]-2*vyf[i+1,j+1]+vyf[i+2,j+1])
            end      
        else()
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf[i-1,j]-2*vyf[i,j]+vyf[i+1,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf[i-1,j+1]-2*vyf[i,j+1]+vyf[i+1,j+1])
            end
        end
        # Compute vy
        vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34
        
        # Change coordinates to obtain B;C;D points
        if(rk==1 || rk==2)
            xcur=xA-dtm/2*vxm[rk]
            ycur=yA-dtm/2*vym[rk]
        elseif[rk==3]
            xcur=xA-dtm*vxm[rk]
            ycur=yA-dtm*vym[rk]
        end
    end

    # Compute effective velocity
    vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
    vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
    # Trace the node backward
    xcur=xA-dtm*vxmeff
    ycur=yA-dtm*vymeff
    # Interpolate nodal property
    # SIGMA'xx; P
    # Define i;j indexes for the upper left node
    j=trunc((xcur-xp[1])/dx)+1
    i=trunc((ycur-yp[1])/dy)+1
    if(j<1)
        j=1
    elseif[j>Nx]
        j=Nx
    end
    if(i<1)
        i=1
    elseif[i>Ny]
        i=Ny
    end
    # Compute distances
    dxmj=xcur-xp[j]
    dymi=ycur-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute nodal total pressure
    pf0[ii,jj]=pf[i,j]*wtmij+pf[i+1,j]*wtmi1j+ pf[i,j+1]*wtmij1+pf[i+1,j+1]*wtmi1j1
end
end


# Add markers to empty areas
marknumold=marknum
mdis=1e30*ones(Nym,Nxm)
mnum=zeros(Nym,Nxm)
mtyp=zeros(Nym,Nxm)
mpor=zeros(Nym,Nxm)
xxm=dxm/2:dxm:xsize-dxm/2
yym=dym/2:dym:ysize-dym/2
for m=1:1:marknum
    
    # Check markers with the nearest nodes
    # Define i;j indexes for the upper left node
    j=trunc((xm[m]-xxm[1])/dxm)+1
    i=trunc((ym[m]-yym[1])/dym)+1
    if(j<1)
        j=1
    elseif[j>Nxm-1]
        j=Nxm-1
    end
    if(i<1)
        i=1
    elseif[i>Nym-1]
        i=Nym-1
    end
    
    # Check nodes
    # i;j Node
    # Compute distance
    dxmj=xm[m]-xxm[j]
    dymi=ym[m]-yym[i]
    dismij=(dxmj^2+dymi^2)^0.5
    if(dismij<mdis[i,j])
        mdis[i,j]=dismij
        mnum[i,j]=m
        mtyp[i,j]=tm[m]
        mpor[i,j]=phim[m]
    end
    # i+1;j Node
    # Compute distance
    dxmj=xm[m]-xxm[j]
    dymi=ym[m]-yym[i+1]
    dismi1j=(dxmj^2+dymi^2)^0.5
    if(dismi1j<mdis[i+1,j])
        mdis[i+1,j]=dismi1j
        mnum[i+1,j]=m
        mtyp[i+1,j]=tm[m]
        mpor[i+1,j]=phim[m]
    end
    # i;j+1 Node
    # Compute distance
    dxmj=xm[m]-xxm[j+1]
    dymi=ym[m]-yym[i]
    dismij1=(dxmj^2+dymi^2)^0.5
    if(dismij1<mdis[i,j+1])
        mdis[i,j+1]=dismij1
        mnum[i,j+1]=m
        mtyp[i,j+1]=tm[m]
        mpor[i,j+1]=phim[m]
    end
    # i+1;j+1 Node
    # Compute distance
    dxmj=xm[m]-xxm[j+1]
    dymi=ym[m]-yym[i+1]
    dismi1j1=(dxmj^2+dymi^2)^0.5
    if(dismi1j1<mdis[i+1,j+1])
        mdis[i+1,j+1]=dismi1j1
        mnum[i+1,j+1]=m
        mtyp[i+1,j+1]=tm[m]
        mpor[i+1,j+1]=phim[m]
    end
end

dii=5*Nxmc
djj=5*Nymc

for j=1:1:Nxm
    for i=1:1:Nym
        if(mnum[i,j]==0)
            # Serch surrounding nodes
            for jj=j-djj:1:j+djj
                for ii=i-dii:1:i+dii
                    if(ii>=1 && ii<=Nym && jj>=1 && jj<=Nxm && mnum[ii,jj]>0)
                        # Compute distance
                        m=mnum[ii,jj]
                        dxmj=xm[m]-xxm[j]
                        dymi=ym[m]-yym[i]
                        dismij=(dxmj^2+dymi^2)^0.5
                        if(dismij<mdis[i,j])
                            mdis[i,j]=dismij
                            mnum[i,j]=-m
                            mtyp[i,j]=-tm[m]
                            mpor[i,j]=phim[m]
                        end
                    end
                end
            end
            # Add New marker
            if(mnum[i,j]<0)
                # Add marker number
                marknum=marknum+1
                # Assign marker coordinates
                xm[marknum]=xxm[j]+(rand()-0.5)*dxm
                ym[marknum]=yym[i]+(rand()-0.5)*dym
                # Copy marker properties
                m=-mnum[i,j]
                tm[marknum]=tm[m]; # Material type()
                tkm[marknum]=tkm[m]; # Temperature
                phim[marknum]=phim[m]; # Porosity
                sxxm[marknum]=sxxm[m]; # SIGMA'xx, Pa
                sxym[marknum]=sxym[m]; # SIGMAxy, Pa
                etavpm[marknum]=etavpm[m]; # Visco-plastic viscosity, Pa
            end
        end
    end
end           
marknumnew=marknum



# Update timesum
global timesum=timesum+dtm

# Translate vx;vy & qxD;qyD into polar coordinates for visualization
vrp=sqrt((vxp+vxp[(Ny+1)*(Nx+1)/2]).^2+(vyp+vyp[(Ny+1)*(Nx+1)/2]).^2)
qrD=sqrt((qxD+qxD[(Ny+1)*(Nx+1)/2]).^2+(qyD+qyD[(Ny+1)*(Nx+1)/2]).^2)

if(trunc(timestep/savematstep)*savematstep==timestep)
    namemat    =  [nname,num2str(timestep)]
    save(namemat)
    fdata=fopen("file.txt','wt")
    fprintf(fdata,"#d",timestep)
    fclose(fdata)
end

# if timesum > 15*3600*24*365.25*1000000
#     break
# end
# end

