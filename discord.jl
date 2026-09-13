using HTTP
using JSON3
using Dates
using Plots

const WEBHOOK = ""


function send_message(msg)
    try
        HTTP.post(
            WEBHOOK,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict("content" => msg))
        )
    catch err
        @warn "Failed to send Discord message" exception=(err, catch_backtrace())
    end
end

function send_image(path; message="")
    try
        open(path) do io
            body = HTTP.Form(Dict(
                "payload_json" => JSON3.write(Dict(
                    "content" => message
                )),
                "file1" => HTTP.Multipart(basename(path), io)
            ))

            HTTP.post(WEBHOOK, [], body)
        end
    catch err
        @warn "Failed to send Discord image" exception=(err, catch_backtrace())
    end
end

function notify_start(label, problem, strategy)

    send_message("""
**Training Started [$label]**

Problem: $problem

Strategy: $(nameof(typeof(strategy)))

Started: $(Dates.now())
""")

end

function notify_progress(label, iter, loss, elapsed)
    send_message("""
**Training Update [$label]**

Iteration: $iter

Loss: $(round(loss, sigdigits = 5))

Elapsed: $(round(elapsed, digits = 1)) s
""")
end

function notify_checkpoint(label, iter)

    send_message("Checkpoint saved at iteration **$iter** **[$label]**")

end

function notify_finish(label, loss, elapsed)

    send_message("""
**Training Finished [$label]**

Final Loss: $(round(loss, sigdigits = 5))

Elapsed: $(round(elapsed, digits = 1)) s
""")

end

function notify_error(label, err)

    send_message("""
**Training Crashed [$label]**

$(typeof(err))

$(err)
""")

end

function notify_loss_plot(label, iter, loss, losses)

    p = plot_loss(losses)

    tmp = tempname() * ".png"

    savefig(p, tmp)

    send_image(
        tmp;
        message = """
**Loss Curve [$label]**

Iteration: $iter

Loss: $(round(loss, sigdigits = 5))
"""
    )

    rm(tmp; force = true)

end

function notify_solution_plot(label, iter, loss, plot)

    tmp = tempname() * ".png"

    savefig(plot, tmp)

    send_image(
        tmp;
        message = """
**Solution [$label]**

Iteration: $iter

Loss: $(round(loss, sigdigits = 5))
"""
    )

    rm(tmp; force = true)

end

function notify_results(label;
    solution = nothing,
    density = nothing,
    loss = nothing,
    residual = nothing
)

    if loss !== nothing
        send_image(loss; message = "**Final Loss Curve [$label]**")
    end

    if solution !== nothing
        send_image(solution; message = "**Final Temperature Solution [$label]**")
    end

    if density !== nothing
        send_image(density; message = "**Final CO Density [$label]**")
    end

    if residual !== nothing
        send_image(residual; message = "**Residual Heatmap [$label]**")
    end

end
