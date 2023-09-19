FROM debian:buster-slim

WORKDIR /bin

COPY ./target/release/umisamedit .

CMD ["/bin/umisamedit"]
