// Built with the production driver globals by potential_table_lookup.py.
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    field_numerical = {4, 3, 2, 1};
    potential_numerical = {16, 9, 4, 1};
    potential_derivative_numerical = {8, 6, 4, 2};
    rescale_B = 1.0;

    const std::string test = argv[1];
    double v, vp;
    if (test == "below" || test == "above") {
        const double field = test == "below" ? std::nextafter(1.0, 0.0)
                                             : std::nextafter(4.0, 5.0);
        evaluate_potential_from_value(field, 3, 1, nullptr, v, vp);
        return 0; // The runner expects the lookup itself to reject this input.
    }

    if (test != "knots") return 2;
    for (double field : field_numerical) {
        for (int hint : {1, 2, 3}) {
            evaluate_potential_from_value(field, hint, 1, nullptr, v, vp);
            if (v != field * field || vp != 2 * field) {
                std::fprintf(stderr, "field=%g hint=%d: V=%g V'=%g\n", field, hint, v, vp);
                return 1;
            }
        }
    }
    return 0;
}
