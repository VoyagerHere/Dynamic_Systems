Threads.@threads for k in eachindex(G_LIST)
  println("Thread $(Threads.threadid()) starting iteration $k")
  # Your code here
  println("Thread $(Threads.threadid()) finished iteration $k")
end