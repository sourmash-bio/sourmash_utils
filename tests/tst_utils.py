import os.path


def get_test_data(filename):
    filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                            filename)
    return filepath

