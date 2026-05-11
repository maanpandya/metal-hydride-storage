using Ferrite, FerriteGmsh, Tensors, LinearAlgebra, SparseArrays, WriteVTK
using FerriteGmsh: Gmsh

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
const nu = 1.0 / 1000.0                 # kinematic viscosity [m^2/s]
const rho = 2.0                         # density [kg/m^3]
const p0_inlet = 1.0                    # prescribed total pressure at inlet [Pa]
const w_oseen = Vec((0.1, 0.0))         # Oseen reference velocity [m/s]

# Total-pressure penalty coefficient
const alpha_bc = 100.0

# Geometry constants (for corner exclusion on inlet)
const H_channel = 0.41
const tol_corner = 0.0011

# Newton settings
const newton_tol = 1e-8
const newton_maxiter = 30

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
function setup_grid(h = 0.05)
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)
    dim = 2

    rect_tag = gmsh.model.occ.add_rectangle(0, 0, 0, 1.1, 0.41)
    circle_tag = gmsh.model.occ.add_circle(0.2, 0.2, 0, 0.05)
    circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
    circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
    gmsh.model.occ.cut([(dim, rect_tag)], [(dim, circle_surf_tag)])
    gmsh.model.occ.synchronize()

    gmsh.model.model.add_physical_group(dim - 1, [6], -1, "bottom")
    gmsh.model.model.add_physical_group(dim - 1, [7], -1, "left")
    gmsh.model.model.add_physical_group(dim - 1, [8], -1, "right")
    gmsh.model.model.add_physical_group(dim - 1, [9], -1, "top")
    gmsh.model.model.add_physical_group(dim - 1, [5], -1, "hole")

    gmsh.option.setNumber("Mesh.Algorithm", 11)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)

    gmsh.model.mesh.generate(dim)
    grid = togrid()
    Gmsh.finalize()
    return grid
end

# ---------------------------------------------------------------------------
# FE spaces and values
# ---------------------------------------------------------------------------
function setup_fevalues(ipv, ipp, ipg)
    qr = QuadratureRule{RefQuadrilateral}(4)
    cvv = CellValues(qr, ipv, ipg)
    cvp = CellValues(qr, ipp, ipg)

    qr_fac = FacetQuadratureRule{RefQuadrilateral}(4)
    fvv = FacetValues(qr_fac, ipv, ipg)
    fvp = FacetValues(qr_fac, ipp, ipg)
    return cvv, cvp, fvv, fvp
end

function setup_dofs(grid, ipv, ipp)
    dh = DofHandler(grid)
    add!(dh, :v, ipv)
    add!(dh, :p, ipp)
    close!(dh)
    return dh
end

# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------
function setup_constraints(dh)
    ch = ConstraintHandler(dh)

    noslip_set = union(
        getfacetset(dh.grid, "top"),
        getfacetset(dh.grid, "bottom"),
        getfacetset(dh.grid, "hole"),
    )
    add!(ch, Dirichlet(:v, noslip_set, (x, t) -> Vec((0.0, 0.0)), [1, 2]))

    # Pressure gauge at outlet to remove null-space.
    add!(ch, Dirichlet(:p, getfacetset(dh.grid, "right"), (x, t) -> 0.0))

    close!(ch)
    update!(ch, 0.0)
    return ch
end

# ---------------------------------------------------------------------------
# Oseen linear stiffness assembly
# ---------------------------------------------------------------------------
function assemble_oseen_linear!(K, dh, cvv, cvp, nu, w)
    assembler = start_assemble(K)
    ke = zeros(ndofs_per_cell(dh), ndofs_per_cell(dh))

    range_v = dof_range(dh, :v)
    ndofs_v = length(range_v)
    range_p = dof_range(dh, :p)
    ndofs_p = length(range_p)

    phi_v = Vector{Vec{2, Float64}}(undef, ndofs_v)
    grad_phi_v = Vector{Tensor{2, 2, Float64, 4}}(undef, ndofs_v)
    div_phi_v = Vector{Float64}(undef, ndofs_v)
    phi_p = Vector{Float64}(undef, ndofs_p)

    for cell in CellIterator(dh)
        Ferrite.reinit!(cvv, cell)
        Ferrite.reinit!(cvp, cell)
        ke .= 0.0

        for qp in 1:getnquadpoints(cvv)
            dOmega = getdetJdV(cvv, qp)

            for i in 1:ndofs_v
                phi_v[i] = shape_value(cvv, qp, i)
                grad_phi_v[i] = shape_gradient(cvv, qp, i)
                div_phi_v[i] = shape_divergence(cvv, qp, i)
            end
            for i in 1:ndofs_p
                phi_p[i] = shape_value(cvp, qp, i)
            end

            for (i, I) in pairs(range_v), (j, J) in pairs(range_v)
                ke[I, J] += (nu * (grad_phi_v[i] ⊡ grad_phi_v[j]) + phi_v[i] ⋅ (grad_phi_v[j] ⋅ w)) * dOmega
            end

            for (i, I) in pairs(range_v), (j, J) in pairs(range_p)
                ke[I, J] += (-div_phi_v[i] * phi_p[j]) * dOmega
            end

            for (i, I) in pairs(range_p), (j, J) in pairs(range_v)
                ke[I, J] += (-phi_p[i] * div_phi_v[j]) * dOmega
            end
        end

        assemble!(assembler, celldofs(cell), ke)
    end

    return K
end

# ---------------------------------------------------------------------------
# Nonlinear inlet total-pressure residual and Jacobian terms
# ---------------------------------------------------------------------------
@inline total_p_residual(ux, p_static) = ux^2 - 2.0 * (p0_inlet - p_static) / rho
@inline dresidual_dux(ux) = 2.0 * ux
@inline dresidual_dp() = 2.0 / rho

function add_inlet_total_pressure_residual!(R, u, dh, fvv, fvp, inlet_boundary)
    v_range = dof_range(dh, :v)
    p_range = dof_range(dh, :p)

    for facet in FacetIterator(dh, inlet_boundary)
        Ferrite.reinit!(fvv, facet)
        Ferrite.reinit!(fvp, facet)

        celld = celldofs(facet)
        v_dofs = @view celld[v_range]
        p_dofs = @view celld[p_range]
        v_e = u[v_dofs]
        p_e = u[p_dofs]

        nphi_v = getnbasefunctions(fvv)
        Re = zeros(length(v_dofs))

        for qp in 1:getnquadpoints(fvv)
            x = spatial_coordinate(fvv, qp, getcoordinates(facet))
            if x[2] <= tol_corner || x[2] >= H_channel - tol_corner
                continue
            end

            dGamma = getdetJdV(fvv, qp)
            ux = function_value(fvv, qp, v_e)[1]
            p_static = function_value(fvp, qp, p_e)
            res = total_p_residual(ux, p_static)

            for i in 1:nphi_v
                phi_x_i = shape_value(fvv, qp, i)[1]
                Re[i] -= alpha_bc * res * phi_x_i * dGamma
            end
        end

        assemble!(R, v_dofs, Re)
    end
end

function add_inlet_total_pressure_jacobian!(J, u, dh, fvv, fvp, inlet_boundary)
    v_range = dof_range(dh, :v)
    p_range = dof_range(dh, :p)

    assembler = start_assemble(J; fillzero = false)

    for facet in FacetIterator(dh, inlet_boundary)
        Ferrite.reinit!(fvv, facet)
        Ferrite.reinit!(fvp, facet)

        celld = celldofs(facet)
        v_dofs = @view celld[v_range]
        p_dofs = @view celld[p_range]
        v_e = u[v_dofs]
        p_e = u[p_dofs]

        nphi_v = getnbasefunctions(fvv)
        nphi_p = getnbasefunctions(fvp)
        Jvv = zeros(length(v_dofs), length(v_dofs))

        for qp in 1:getnquadpoints(fvv)
            x = spatial_coordinate(fvv, qp, getcoordinates(facet))
            if x[2] <= tol_corner || x[2] >= H_channel - tol_corner
                continue
            end

            dGamma = getdetJdV(fvv, qp)
            ux = function_value(fvv, qp, v_e)[1]
            dres_dux = dresidual_dux(ux)
            dres_dp = dresidual_dp()

            for i in 1:nphi_v
                phi_x_i = shape_value(fvv, qp, i)[1]

                for j in 1:nphi_v
                    phi_x_j = shape_value(fvv, qp, j)[1]
                    Jvv[i, j] -= alpha_bc * dres_dux * phi_x_i * phi_x_j * dGamma
                end

                global_i = v_dofs[i]
                for j in 1:nphi_p
                    psi_j = shape_value(fvp, qp, j)
                    global_j = p_dofs[j]
                    J[global_i, global_j] -= alpha_bc * dres_dp * phi_x_i * psi_j * dGamma
                end
            end
        end

        assemble!(assembler, v_dofs, Jvv)
    end
end

# ---------------------------------------------------------------------------
# Nonlinear solve: F(u)=0 with Newton
# ---------------------------------------------------------------------------
function solve_oseen_total_pressure(dh, ch, K_lin, fvv, fvp, inlet_boundary)
    u = zeros(ndofs(dh))
    apply!(u, ch)

    R = zeros(ndofs(dh))

    for it in 1:newton_maxiter
        # Residual
        mul!(R, K_lin, u)
        add_inlet_total_pressure_residual!(R, u, dh, fvv, fvp, inlet_boundary)

        # Jacobian
        J = copy(K_lin)
        add_inlet_total_pressure_jacobian!(J, u, dh, fvv, fvp, inlet_boundary)

        # Enforce Dirichlet constraints in the Newton system
        apply!(J, R, ch)

        res_norm = norm(R)
        println("Newton iteration $it: ||R|| = $(res_norm)")
        if res_norm < newton_tol
            println("Newton converged.")
            return u
        end

        du = -(J \ R)
        u .+= du
        apply!(u, ch)
    end

    error("Newton did not converge in $(newton_maxiter) iterations.")
end

# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
function main()
    println("Setting up grid...")
    grid = setup_grid(0.05)

    ipv = Lagrange{RefQuadrilateral, 2}()^2
    ipp = Lagrange{RefQuadrilateral, 1}()
    ipg = Lagrange{RefQuadrilateral, 1}()

    dh = setup_dofs(grid, ipv, ipp)
    cvv, cvp, fvv, fvp = setup_fevalues(ipv, ipp, ipg)
    ch = setup_constraints(dh)

    inlet_boundary = getfacetset(grid, "left")

    println("Assembling Oseen linear matrix (nu = $nu, w = $w_oseen)...")
    coupling = [true true; true false]
    K_lin = allocate_matrix(dh, ch; coupling = coupling)
    assemble_oseen_linear!(K_lin, dh, cvv, cvp, nu, w_oseen)

    println("Solving nonlinear Oseen system with inlet total-pressure BC...")
    u = solve_oseen_total_pressure(dh, ch, K_lin, fvv, fvp, inlet_boundary)

    println("Writing VTK output...")
    VTKGridFile("oseen-p0-inlet-cyl", grid) do vtk
        write_solution(vtk, dh, u)
    end

    println("Done! Solution written to oseen-p0-inlet-cyl.vtu")
end

main()
