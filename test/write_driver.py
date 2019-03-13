import FRUIT


def main():
    test_modules = ["linalg_test.f90", "l2m_test.f90"]
    driver = "test_driver.f90"
    suite = FRUIT.test_suite(test_modules)
    suite.write(driver)
    return


if __name__ == '__main__':
    main()
