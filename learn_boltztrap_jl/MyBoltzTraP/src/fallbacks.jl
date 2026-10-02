# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Fallback errors

function savefig_publication(args...; kwargs...)
    throw(
        ArgumentError(
            "savefig_publication requires PythonPlot.jl. Run `using PythonPlot` first.",
        ),
    )
end
