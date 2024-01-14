#include "../formatter_fmt.h"
#include <gtest/gtest.h>

/************************************************
 *  unit test of class Fmt
 ***********************************************/

/**
 * - Tested Functions:
 *   - Fmt()
 *     - default constructor, default values see formatter.h
 *   - Fmt(int width, int precision, char fillChar, bool fixed, bool right, bool error)
 *     - parameterized constructor
 *   - set_width(int width)
 *     - setter of width
 *   - set_precision(int precision)
 *     - setter of precision
 *   - set_fillChar(char fillChar)
 *     - setter of fillChar
 *   - set_fixed(bool fixed)
 *     - setter of scientific notation boolean, fixed means not scientific
 *   - set_right(bool right)
 *     - setter of right alignment boolean
 *   - set_error(bool error)
 *     - setter of + sign for positive numbers boolean
 *   - format(T value)
 *     - format a value, according to width, precision, fillChar, fixed, right, error
 */

class FmtTest : public testing::Test
{
};

TEST_F(FmtTest, DefaultConstructor) {
    formatter::Fmt fmt;
    EXPECT_EQ(fmt.get_width(), 4);
    EXPECT_EQ(fmt.get_precision(), 2);
    EXPECT_EQ(fmt.get_fillChar(), ' ');
    EXPECT_EQ(fmt.get_fixed(), true);
    EXPECT_EQ(fmt.get_right(), true);
    EXPECT_EQ(fmt.get_error(), false);
}

TEST_F(FmtTest, ParameterizedConstructor) {
    formatter::Fmt fmt(10, 5, '0', false, false, true);
    EXPECT_EQ(fmt.get_width(), 10);
    EXPECT_EQ(fmt.get_precision(), 5);
    EXPECT_EQ(fmt.get_fillChar(), '0');
    EXPECT_EQ(fmt.get_fixed(), false);
    EXPECT_EQ(fmt.get_right(), false);
    EXPECT_EQ(fmt.get_error(), true);
}

TEST_F(FmtTest, Setters) {
    formatter::Fmt fmt;
    fmt.set_width(10);
    fmt.set_precision(5);
    fmt.set_fillChar('0');
    fmt.set_fixed(false);
    fmt.set_right(false);
    fmt.set_error(true);
    EXPECT_EQ(fmt.get_width(), 10);
    EXPECT_EQ(fmt.get_precision(), 5);
    EXPECT_EQ(fmt.get_fillChar(), '0');
    EXPECT_EQ(fmt.get_fixed(), false);
    EXPECT_EQ(fmt.get_right(), false);
    EXPECT_EQ(fmt.get_error(), true);
}

TEST_F(FmtTest, Format) {
    formatter::Fmt fmt;
    EXPECT_EQ(fmt.format(1), "1.00");
    EXPECT_EQ(fmt.format(1.23456789), "1.23");
    EXPECT_EQ(fmt.format(1.23456789e-10), "0.00");
    EXPECT_EQ(fmt.format(1.23456789e10), "12345678900.00");
    EXPECT_EQ(fmt.format(-1), "-1.00");
    EXPECT_EQ(fmt.format(-1.23456789), "-1.23");
    EXPECT_EQ(fmt.format(-1.23456789e-10), "-0.00");
    EXPECT_EQ(fmt.format(-1.23456789e10), "-12345678900.00");

    EXPECT_EQ(fmt.format((std::string)"hello"), "hello");
}

TEST_F(FmtTest, SpecialFormats) {
    formatter::Fmt fmt;
    fmt.set_width(8); fmt.set_precision(4); fmt.set_fillChar(' '); fmt.set_fixed(true); fmt.set_right(true); fmt.set_error(false);
    EXPECT_EQ(fmt.format(1), "  1.0000");
    fmt.set_error(true);
    EXPECT_EQ(fmt.format(1), " +1.0000");
    EXPECT_EQ(fmt.format(-1), " -1.0000");
}

TEST_F(FmtTest, ComplexFormat) {
    formatter::Fmt fmt;
    fmt.set_width(8); fmt.set_precision(4); fmt.set_fillChar(' '); fmt.set_fixed(true); fmt.set_right(true); fmt.set_error(false);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "(  1.0000,  2.0000)");
    fmt.set_error(true);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "( +1.0000, +2.0000)");
    fmt.set_width(20); fmt.set_precision(10); fmt.set_fillChar(' '); fmt.set_fixed(true); fmt.set_right(true); fmt.set_error(false);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "(        1.0000000000,        2.0000000000)");
    fmt.set_error(true);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "(       +1.0000000000,       +2.0000000000)");
    EXPECT_EQ(fmt.format(std::complex<double>(1, -2)), "(       +1.0000000000,       -2.0000000000)");

}