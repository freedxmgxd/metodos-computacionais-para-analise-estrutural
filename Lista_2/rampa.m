function [fy_1, m_1, fy_2, m_2] = rampa(L, w_max, w_min)
    % Ramp part
    fy_1_ramp = (- 7 * (w_max - w_min) * L) / 20;
    m_1_ramp = (- (w_max - w_min) * L * L) / 20;
    fy_2_ramp = (- 3 * (w_max - w_min) * L) / 20;
    m_2_ramp = ((w_max - w_min) * L * L) / 30;

    % Rect part
    fy_1_rect = (- w_min * L) / 2;
    m_1_rect = (- w_min * L * L) / 12;
    fy_2_rect = (- w_min * L) / 2;
    m_2_rect = (w_min * L * L) / 12;

    % Total
    fy_1 = fy_1_ramp + fy_1_rect;
    m_1 = m_1_ramp + m_1_rect;
    fy_2 = fy_2_ramp + fy_2_rect;
    m_2 = m_2_ramp + m_2_rect;
end

