#include "socket_driver.h"

#include "source_relax/socket_driver_handlers.h"
#include "source_relax/socket_driver_utils.h"
#include "source_base/timer.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver.h"
#include "source_io/module_parameter/input_parameter.h"

#include <exception>
#include <fstream>
#include <string>

using SocketDriverUtils::ComputedFrame;
using SocketDriverUtils::DriverState;
using SocketDriverUtils::ipi_cell_bohr_from_unitcell;
using SocketDriverUtils::is_root;
using SocketDriverUtils::quit_if_root_io_failed;
using SocketDriverUtils::socket_address;
using SocketDriverHandlers::DriverContext;

namespace
{
void connect_on_root(IpiSocket& socket, std::ofstream& ofs_running)
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
}

void log_peer_closed(std::ofstream& ofs_running)
{
    if (is_root())
    {
        ofs_running << " ABACUS socket driver exiting after peer closed connection" << std::endl;
    }
}
} // namespace

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
        connect_on_root(socket, ofs_running);

        DriverContext context;
        context.esolver = p_esolver;
        context.ucell = &ucell;
        context.inp = &inp;
        context.state = DriverState::NeedInit;
        context.nat_return = ucell.nat;
        context.reference_cell = ipi_cell_bohr_from_unitcell(ucell);

        while (true)
        {
            const std::string header
                = SocketDriverHandlers::read_header_bcast(socket, context.state);
            if (header.empty())
            {
                log_peer_closed(ofs_running);
                break;
            }
            else if (header == "STATUS")
            {
                SocketDriverHandlers::handle_status(socket, context.state);
            }
            else if (header == "INIT")
            {
                SocketDriverHandlers::handle_init(socket, context, ofs_running);
            }
            else if (header == "POSDATA")
            {
                SocketDriverHandlers::handle_posdata(socket, context, ofs_running);
            }
            else if (header == "GETFORCE")
            {
                SocketDriverHandlers::handle_getforce(socket, context);
            }
            else if (header == "EXIT")
            {
                SocketDriverHandlers::handle_exit(ofs_running);
                break;
            }
            else
            {
                quit_if_root_io_failed(is_root() ? 1 : 0,
                                       is_root() ? "unknown i-PI header: " + header : "");
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
