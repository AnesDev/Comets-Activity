using HTTP
using JSON3
using Dates
using Plots

const WEBHOOK = "https://discord.com/api/webhooks/1524517919369859163/k7AT4QY-UK6CuNG8LgXPxwDdZJo5dvWLhhA2Lq-wjlzOTXb5orA65KXWIh28HQwDjR8U"


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
function notify_start(problem, strategy)

    send_message("""
**Training Started**

Problem: $problem

Strategy: $(nameof(typeof(strategy)))

Started: $(Dates.now())
""")

end

function notify_progress(iter, loss, elapsed)

    send_message("""
**Training Update**

Iteration: $iter

Loss: $(round(loss, sigdigits = 5))

Elapsed: $(round(elapsed, digits = 1)) s
""")

end

function notify_checkpoint(iter)

    send_message("Checkpoint saved at iteration **$iter**")

end

function notify_finish(loss, elapsed)

    send_message("""
**Training Finished**

Final Loss: $(round(loss, sigdigits = 5))

Elapsed: $(round(elapsed, digits = 1)) s
""")

end

function notify_error(err)

    send_message("""
**Training Crashed**

$(typeof(err))

$(err)
""")

end

function notify_loss_plot(iter, loss, losses)

    p = plot_loss(losses)

    tmp = tempname() * ".png"

    savefig(p, tmp)

    send_image(
        tmp;
        message = """
**Loss Curve**

Iteration: $iter

Loss: $(round(loss, sigdigits = 5))
"""
    )

    rm(tmp; force = true)

end

function notify_solution_plot(iter, loss, plot)

    tmp = tempname() * ".png"

    savefig(plot, tmp)

    send_image(
        tmp;
        message = """
**Solution**

Iteration: $iter

Loss: $(round(loss, sigdigits = 5))
"""
    )

    rm(tmp; force = true)

end
function notify_results(;
    solution = nothing,
    density = nothing,
    loss = nothing,
    residual = nothing
)

    if loss !== nothing
        send_image(loss; message = "**Final Loss Curve**")
    end

    if solution !== nothing
        send_image(solution; message = "**Final Temperature Solution**")
    end

    if density !== nothing
        send_image(density; message = "**Final CO Density**")
    end

    if residual !== nothing
        send_image(residual; message = "**Residual Heatmap**")
    end

end

