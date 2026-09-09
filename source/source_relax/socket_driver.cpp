#include "socket_driver.h"

#include "source_relax/socket_driver_utils.h"
#include "source_relax/socket_frame.h"
#include "source_relax/socket_ipi.h"
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace SocketDriverUtils;

void Socket_Driver::socket_driver(ModuleESolver::ESolver* p_esolver,
                                  UnitCell& ucell,
                                  const Input_para& inp,
                                  std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Socket_Driver", "socket_driver");
    ModuleBase::timer::start("Socket_Driver", "socket_driver");

    if (p_esolver == nullptr)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", "socket driver requires a valid ESolver.");
    }
    IpiSocket socket;

    try
    {
        int io_failed = 0;
        std::string io_message;
        if (is_root())
        {
            try
            {
                const std::string address = socket_address();
                ofs_running << " ABACUS socket driver connecting to i-PI endpoint " << address << std::endl;
                socket.connect(address);
            }
            catch (const std::exception& exc)
            {
                io_failed = 1;
                io_message = exc.what();
            }
        }
        quit_if_root_io_failed(io_failed, io_message);

        DriverState state = DriverState::NeedInit;
        int istep = 0;
        const int nat_return = ucell.nat;
        ComputedFrame published;

        const std::vector<double> reference_cell = ipi_cell_bohr_from_unitcell(ucell);
        bool checked_initial_positions = false;

        while (true)
        {
            std::string header;
            io_failed = 0;
            io_message.clear();
            if (is_root())
            {
                try
                {
                    header = socket.read_header();
                }
                catch (const IpiSocketClosed&)
                {
                    if (state == DriverState::HasData)
                    {
                        io_failed = 1;
                        io_message = "i-PI peer closed while a computed frame was pending";
                    }
                    else
                    {
                        header.clear();
                    }
                }
                catch (const std::exception& exc)
                {
                    io_failed = 1;
                    io_message = exc.what();
                }
            }
            quit_if_root_io_failed(io_failed, io_message);
            header = bcast_header(header);

            if (header.empty())
            {
                if (is_root())
                {
                    ofs_running << " ABACUS socket driver exiting after peer closed connection" << std::endl;
                }
                break;
            }
            else if (header == "STATUS")
            {
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    try
                    {
                        if (state == DriverState::HasData)
                        {
                            socket.write_header("HAVEDATA");
                        }
                        else if (state == DriverState::Ready)
                        {
                            socket.write_header("READY");
                        }
                        else
                        {
                            socket.write_header("NEEDINIT");
                        }
                    }
                    catch (const std::exception& exc)
                    {
                        io_failed = 1;
                        io_message = exc.what();
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
            }
            else if (header == "INIT")
            {
                std::int32_t rid = 0;
                std::int32_t nbytes = 0;
                std::string params;
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    if (state != DriverState::NeedInit)
                    {
                        io_failed = 1;
                        io_message = "INIT requires NEEDINIT state";
                    }
                    else
                    {
                        try
                        {
                            rid = socket.read_int32();
                            nbytes = socket.read_int32();
                            if (nbytes < 0)
                            {
                                io_failed = 1;
                                io_message = "negative INIT payload length from i-PI socket";
                            }
                            else if (nbytes > MAX_INIT_BYTES)
                            {
                                io_failed = 1;
                                io_message = "INIT payload exceeds the 1 MiB socket limit";
                            }
                            else if (nbytes > 0)
                            {
                                params = socket.read_string(static_cast<std::size_t>(nbytes));
                            }
                        }
                        catch (const std::exception& exc)
                        {
                            io_failed = 1;
                            io_message = exc.what();
                        }
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
                bcast_socket_int32(rid);
                bcast_socket_int32(nbytes);
                if (nbytes > 0 && is_root())
                {
                    ofs_running << " ABACUS socket INIT params bytes " << nbytes << std::endl;
                }
                state = DriverState::Ready;
                if (is_root())
                {
                    ofs_running << " ABACUS socket INIT replica " << rid << std::endl;
                }
            }
            else if (header == "POSDATA")
            {
                SocketFrame::Matrix9 cell = {{0.0}};
                SocketFrame::Matrix9 inv_cell = {{0.0}};
                std::int32_t nat_socket = 0;
                std::vector<double> positions;
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    if (state != DriverState::Ready)
                    {
                        io_failed = 1;
                        io_message = "POSDATA requires READY state";
                    }
                    else
                    {
                        try
                        {
                            const std::vector<double> cell_values = socket.read_doubles(9);
                            const std::vector<double> inverse_values = socket.read_doubles(9);
                            std::copy(cell_values.begin(), cell_values.end(), cell.begin());
                            std::copy(inverse_values.begin(), inverse_values.end(), inv_cell.begin());
                            nat_socket = socket.read_int32();
                            SocketFrame::CellValidation validation
                                = SocketFrame::validate_ipi_cell(cell,
                                                                 inv_cell,
                                                                 MAX_CELL_CONDITION,
                                                                 INVERSE_ABSOLUTE_TOLERANCE,
                                                                 INVERSE_RELATIVE_TOLERANCE);
                            if (!validation.ok)
                            {
                                io_failed = 1;
                                io_message = "invalid POSDATA cell: " + validation.message;
                            }
                            std::size_t coordinate_count = 0;
                            if (io_failed == 0
                                && !SocketFrame::checked_position_count(nat_socket,
                                                                        ucell.nat,
                                                                        coordinate_count,
                                                                        io_message))
                            {
                                io_failed = 1;
                            }
                            if (io_failed == 0)
                            {
                                positions = socket.read_doubles(coordinate_count);
                                if (!SocketFrame::validate_positions(positions,
                                                                     coordinate_count,
                                                                     io_message))
                                {
                                    io_failed = 1;
                                }
                            }
                        }
                        catch (const std::exception& exc)
                        {
                            io_failed = 1;
                            io_message = exc.what();
                        }
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
                bcast_socket_int32(nat_socket);
                std::vector<double> cell_values(cell.begin(), cell.end());
                std::vector<double> inverse_values(inv_cell.begin(), inv_cell.end());
                bcast_double_vector(cell_values);
                bcast_double_vector(inverse_values);
                if (!is_root())
                {
                    cell = {{0.0}};
                    inv_cell = {{0.0}};
                    std::copy(cell_values.begin(), cell_values.end(), cell.begin());
                    std::copy(inverse_values.begin(), inverse_values.end(), inv_cell.begin());
                    if (nat_socket >= 0)
                    {
                        positions.assign(static_cast<std::size_t>(3 * nat_socket), 0.0);
                    }
                }
                bcast_double_vector(positions);

                const double max_cell_delta_bohr = max_abs_delta(std::vector<double>(cell.begin(), cell.end()), reference_cell);
                if (max_cell_delta_bohr > unchanged_cell_tolerance(cell))
                {
                    ModuleBase::WARNING_QUIT("ABACUS socket", "variable-cell socket updates are not supported yet.");
                }
                if (!checked_initial_positions)
                {
                    checked_initial_positions = true;
                    if (max_wrapped_direct_delta_from_unitcell(ucell, positions) > 1.0e-5 && is_root())
                    {
                        ModuleBase::WARNING(
                            "ABACUS socket",
                            "first POSDATA positions are not PBC-equivalent to STRU atom order; "
                            "i-PI POSDATA carries no species, so the client atoms should use the same atom order as STRU.");
                    }
                }

                try
                {
                    set_positions_from_ipi_bohr(ucell, positions);
                }
                catch (const std::exception& exc)
                {
                    fail_during_collective_stage("set_positions", exc.what());
                }
                catch (...)
                {
                    fail_during_collective_stage("set_positions",
                                                 "unknown socket position update failure");
                }
                try
                {
                    p_esolver->runner(ucell, istep);
                }
                catch (const std::exception& exc)
                {
                    fail_during_collective_stage("runner", exc.what());
                }
                catch (...)
                {
                    fail_during_collective_stage("runner",
                                                 "unknown socket runner failure");
                }
                ComputedFrame computed;
                computed.scf_converged = all_ranks_converged(p_esolver->conv_esolver);
                if (!computed.scf_converged && is_root())
                {
                    ModuleBase::WARNING(
                        "ABACUS socket",
                        "SCF did not converge; returning the available frame and marking it in i-PI extras.");
                }
                double energy_ry = 0.0;
                try
                {
                    energy_ry = p_esolver->cal_energy();
                }
                catch (const std::exception& exc)
                {
                    fail_during_collective_stage("cal_energy", exc.what());
                }
                catch (...)
                {
                    fail_during_collective_stage("cal_energy",
                                                 "unknown socket energy failure");
                }
                int local_failed = std::isfinite(energy_ry) ? 0 : 1;
                throw_if_any_rank_failed(local_failed,
                                         local_failed == 0 ? "" : "socket energy is not finite");
                if (!std::isfinite(energy_ry))
                {
                    ModuleBase::WARNING_QUIT("ABACUS socket", "socket energy is not finite.");
                }
                computed.energy_hartree = energy_ry * RY_TO_HARTREE;
                if (is_root())
                {
                    ofs_running << " ABACUS socket return energy "
                                << energy_ry << " Ry, "
                                << energy_ry * ModuleBase::Ry_to_eV << " eV, "
                                << computed.energy_hartree << " Ha" << std::endl;
                }
                ModuleBase::matrix force;
                if (inp.cal_force)
                {
                    try
                    {
                        p_esolver->cal_force(ucell, force);
                    }
                    catch (const std::exception& exc)
                    {
                        fail_during_collective_stage("cal_force", exc.what());
                    }
                    catch (...)
                    {
                        fail_during_collective_stage("cal_force",
                                                     "unknown socket force failure");
                    }
                    local_failed = 0;
                    std::string local_message;
                    try
                    {
                        computed.forces_hartree_per_bohr = flatten_forces_hartree_per_bohr(force, ucell.nat);
                    }
                    catch (const std::exception& exc)
                    {
                        local_failed = 1;
                        local_message = exc.what();
                    }
                    catch (...)
                    {
                        local_failed = 1;
                        local_message = "unknown socket force validation failure";
                    }
                    throw_if_any_rank_failed(local_failed, local_message);
                    computed.forces_present = true;
                }
                if (inp.cal_stress)
                {
                    ModuleBase::matrix stress;
                    try
                    {
                        p_esolver->cal_stress(ucell, stress);
                    }
                    catch (const std::exception& exc)
                    {
                        fail_during_collective_stage("cal_stress", exc.what());
                    }
                    catch (...)
                    {
                        fail_during_collective_stage("cal_stress",
                                                     "unknown socket stress failure");
                    }
                    local_failed = 0;
                    std::string local_message;
                    try
                    {
                        const SocketFrame::VirialConversion virial
                            = SocketFrame::make_ipi_virial(matrix9_from_stress(stress),
                                                           ucell.omega,
                                                           STRESS_ABSOLUTE_TOLERANCE,
                                                           STRESS_RELATIVE_TOLERANCE);
                        if (!virial.ok)
                        {
                            throw std::runtime_error(virial.message);
                        }
                        computed.virial_wire_hartree = virial.wire_virial_hartree;
                    }
                    catch (const std::exception& exc)
                    {
                        local_failed = 1;
                        local_message = exc.what();
                    }
                    catch (...)
                    {
                        local_failed = 1;
                        local_message = "unknown socket stress validation failure";
                    }
                    throw_if_any_rank_failed(local_failed, local_message);
                    computed.stress_present = true;
                }
                computed.valid = true;
                published = computed;
                ++istep;
                state = DriverState::HasData;
            }
            else if (header == "GETFORCE")
            {
                io_failed = 0;
                io_message.clear();
                if (is_root())
                {
                    try
                    {
                        if (state != DriverState::HasData || !published.valid)
                        {
                            throw std::runtime_error("GETFORCE requires HAVEDATA state and a valid frame");
                        }
                        socket.write_header("FORCEREADY");
                        socket.write_double(published.energy_hartree);
                        socket.write_int32(static_cast<std::int32_t>(nat_return));
                        const std::vector<double> forces
                            = published.forces_present
                                  ? published.forces_hartree_per_bohr
                                  : std::vector<double>(static_cast<std::size_t>(3 * nat_return), 0.0);
                        socket.write_doubles(forces);
                        socket.write_doubles(vector_from_matrix9(published.virial_wire_hartree));
                        const std::string extra = properties_extra(published);
                        if (extra.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
                        {
                            throw std::overflow_error("i-PI extras payload is larger than int32");
                        }
                        socket.write_int32(static_cast<std::int32_t>(extra.size()));
                        socket.write_string(extra);
                    }
                    catch (const std::exception& exc)
                    {
                        io_failed = 1;
                        io_message = exc.what();
                    }
                }
                quit_if_root_io_failed(io_failed, io_message);
                published = ComputedFrame();
                state = DriverState::Ready;
            }
            else if (header == "EXIT")
            {
                if (is_root())
                {
                    ofs_running << " ABACUS socket driver received i-PI EXIT" << std::endl;
                }
                break;
            }
            else
            {
                if (is_root())
                {
                    io_failed = 1;
                    io_message = "unknown i-PI header: " + header;
                }
                quit_if_root_io_failed(io_failed, io_message);
            }
        }
    }
    catch (const std::exception& exc)
    {
        ModuleBase::WARNING_QUIT("ABACUS socket", exc.what());
    }

    if (is_root())
    {
        socket.close();
    }

    ModuleBase::timer::end("Socket_Driver", "socket_driver");
}
