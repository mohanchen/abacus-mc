#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/read_input.h"
#include "source_io/module_parameter/system_parameter.h"
/************************************************
 *  unit test of read_input_test_item.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Item_test:
 *     - read in specific values for some items
 */

class InputTest : public testing::Test
{
  protected:
    std::vector<std::pair<std::string, ModuleIO::Input_Item>>::iterator find_label(
        const std::string& label,
        std::vector<std::pair<std::string, ModuleIO::Input_Item>>& bcastfuncs)
    {
        auto it = std::find_if(
            bcastfuncs.begin(),
            bcastfuncs.end(),
            [&label](const std::pair<std::string, ModuleIO::Input_Item>& item) { return item.first == label; });
        return it;
    }
};
ModuleIO::ReadInput readinput(0);
Input_para input;
System_para sys;
std::string output = "";

TEST_F(InputTest, Item_test)
{
    readinput.check_ntype_flag = false;

    {
        input.suffix = "test";
        readinput.set_global_dir(input, sys);

        EXPECT_EQ(sys.global_out_dir, "OUT.test/");
        EXPECT_EQ(sys.global_stru_dir, "OUT.test/STRU/");
        EXPECT_EQ(sys.global_matrix_dir, "OUT.test/matrix/");

        readinput.set_globalv(input, sys);

        input.basis_type = "lcao";
        input.gamma_only = true;
        input.esolver_type = "tddft";
        input.nspin = 2;
        readinput.set_globalv(input, sys);
        EXPECT_EQ(sys.gamma_only_local, 0);

        input.deepks_scf = true;
        input.deepks_out_labels = true;
        readinput.set_globalv(input, sys);
        EXPECT_EQ(sys.deepks_setorb, 1);

        input.nspin = 4;
        input.noncolin = true;
        readinput.set_globalv(input, sys);
        EXPECT_EQ(sys.domag, 1);
        EXPECT_EQ(sys.domag_z, 0);
        EXPECT_EQ(sys.npol, 2);

        input.nspin = 1;
        input.lspinorb = true;
        input.noncolin = false;
        readinput.set_globalv(input, sys);
        EXPECT_EQ(sys.domag, 0);
        EXPECT_EQ(sys.domag_z, 0);
        EXPECT_EQ(sys.npol, 1);
    }
}
