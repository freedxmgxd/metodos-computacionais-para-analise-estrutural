function [fy_1, m_1, fy_2, m_2] = rampa(L, w_max, w_min)
    fy_1 = -7*(w_max - w_min)*L/20 - w_min*L/2;
    m_1 = -(w_max - w_min) *L*L/20 - w_min*L*L/12;
    fy_2 = -3*(w_max - w_min)*L/20 - w_min*L/2;
    m_2 = (w_max - w_min) *L*L/30 + w_min*L*L/12;
end

