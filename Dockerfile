FROM ubuntu

COPY ./target/release/samedit /samedit
RUN chmod +x /samedit

ENTRYPOINT [ "/samedit" ]
