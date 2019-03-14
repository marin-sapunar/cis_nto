#!/usr/bin/python3
""" Generate driver program for unit testing using FRUIT. """
import FRUIT


TEST_MODULES = [
                "linalg_test.f90",
                "assignment_test.f90",
                "l2m_test.f90",
                ]


TEST_DRIVER = "fruit_driver.f90"


def main():
    """ Write the program to TEST_DRIVER file. """
    suite = FRUIT.test_suite(TEST_MODULES)
    suite.write(TEST_DRIVER)
    return


if __name__ == '__main__':
    main()
