

hailctl dataproc start haplo1 \
    -m n1-highmem-4 \
    --worker n1-standard-4 \
    --region us-central1 \
    --max-idle 3h \
    --network dataproc \
    --autoscaling-policy max-50
