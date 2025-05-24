#pragma once
class RunningStat
{
public:
    RunningStat();
    void Clear();
    void Push(double x);
    int NumDataValues();
    double Mean();
    double Variance();
private:
    int m_n;
    double m_oldM, m_newM, m_oldS, m_newS;
};
