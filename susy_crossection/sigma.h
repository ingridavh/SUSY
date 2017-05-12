#ifndef SIGMA_H
#define SIGMA_H


class sigma
{
public:
    sigma(double s);
    double integrate_simpson(double x_max, double x_min, double N);
    double sigma_diff(double x1, double x2, double s_hat);
    double f(double x);
    void same_quarks();
    void different_quarks();
    double m_s;
    double m_mg;
    double m_mq;
    double m_gs;
    double m_m2;
    double m_alpha;
    const double m_pi;
    int m_i;
    int m_error_count;
    int m_no_prob;
};

#endif // SIGMA_H
